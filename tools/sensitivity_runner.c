#define SOLVER_REAL_AS double
#include "config/config.h"
#include "model/input.h"
#include "model/model.h"
#include "simulation/simulation.h"
#include "util/fpclass.h"

#include <framework/batch_runner.h>
#include <framework/sensitivity.h>

#include <math.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* ------------------------------------------------------------------ limits */

#define MAX_PARAMS   20
#define MAX_WORKERS  16

/* ------------------------------------------------------------------ types  */

typedef enum { METRIC_RGR = 0, METRIC_BIOMASS, METRIC_SUCROSE, METRIC_STARCH, METRIC_CARB } MetricKind;
typedef enum { ANALYSIS_MORRIS = 1, ANALYSIS_SOBOL = 2, ANALYSIS_BOTH = 3 } AnalysisKind;

typedef struct {
    Input            base_input;
    SimulationConfig cfg;
    size_t           offsets[MAX_PARAMS];
    MetricKind       metric;
    /* Divergence accounting. The Saltelli estimator relies on row pairing, so a
     * diverged evaluation cannot be deleted from the design; it must be counted
     * and reported alongside the indices. Incremented from worker threads, hence
     * atomic. n_eval counts every call, n_diverged those whose final biomass is
     * non-finite or non-positive. */
    _Atomic unsigned long n_eval;
    _Atomic unsigned long n_diverged;
    _Atomic unsigned long n_bad_endpoint;
} PlantSensCtx;

typedef struct {
    int    photoperiod_h;
    double lit_rgr;
} Scenario;

/* ------------------------------------------------------------------ model evaluate */

/*
 * Receives physically-scaled parameter values from the sensitivity library
 * (sensitivity_evaluate_unit_samples_batch and sensitivity_morris_analyze both
 * scale unit [0,1] samples to physical ranges before calling this function).
 */
static int plant_evaluate(
    const real_t* parameters,
    size_t        parameter_count,
    real_t*       out_value,
    void*         ctx)
{
    PlantSensCtx* pctx = (PlantSensCtx*)ctx;

    Input input = pctx->base_input;
    for (size_t i = 0; i < parameter_count; ++i) {
        real_t* field = (real_t*)((char*)&input + pctx->offsets[i]);
        *field = parameters[i];
    }

    real_t fw0 = input.growth.total_biomass;
    SimulationResult result;
    simulate_days(&pctx->cfg, &input, &result);

    /* Classify before the metric switch, so the count is metric-independent.
     * simulate_days sets result.diverged when the NaN guard tripped at any
     * sampled point OR the step budget was exhausted; that flag, not the final
     * state, is the authoritative classification, because a trajectory can pass
     * through a non-finite state and still land on a finite endpoint that the
     * clamp below would silently accept. Both are counted: n_diverged is the
     * flag, n_bad_endpoint the subset whose endpoint is also unusable. */
    atomic_fetch_add(&pctx->n_eval, 1UL);
    if (result.diverged)
        atomic_fetch_add(&pctx->n_diverged, 1UL);
    /* pm_isnonfinite, not isfinite: this project builds with -ffast-math, under
     * which the library macro folds to a constant and the test disappears. The
     * same trap is documented at include/util/fpclass.h and in the NaN guard. */
    if (pm_isnonfinite((double)input.growth.total_biomass)
        || !(input.growth.total_biomass > REAL_EPSILON))
        atomic_fetch_add(&pctx->n_bad_endpoint, 1UL);

    simulation_result_free(&result);

    switch (pctx->metric) {
        case METRIC_RGR: {
            /* Clamp so log() does not produce NaN/inf for dead plants. */
            real_t fw_end = input.growth.total_biomass;
            if (!(fw_end > REAL_EPSILON))
                fw_end = REAL_EPSILON;
            *out_value = compute_rgr(fw0, fw_end, (real_t)pctx->cfg.days);
            break;
        }
        case METRIC_BIOMASS:
            *out_value = RMAX(input.growth.total_biomass, REAL(0.0));
            break;
        case METRIC_SUCROSE:
            *out_value = RMAX(input.carbohydrates.sucrose, REAL(0.0));
            break;
        case METRIC_STARCH:
            *out_value = RMAX(input.carbohydrates.starch, REAL(0.0));
            break;
        case METRIC_CARB: {
            /* J_carb = 0.5*(sucrose/max_sucrose + starch/max_starch): combined
             * normalised carbohydrate status; range ≈ [0,1]. */
            real_t s_max  = pctx->base_input.carbohydrates.max_sucrose;
            real_t st_max = pctx->base_input.carbohydrates.max_starch;
            real_t s_frac  = RMAX(input.carbohydrates.sucrose, REAL(0.0)) / s_max;
            real_t st_frac = RMAX(input.carbohydrates.starch,  REAL(0.0)) / st_max;
            *out_value = REAL(0.5) * (s_frac + st_frac);
            break;
        }
        default:
            *out_value = REAL(0.0);
            break;
    }
    return 0;
}

/* ==================================================================
 * Deterministic RNG (xorshift64*) — used for Morris trajectory design
 * and for bootstrap resampling.  Self-contained so that every reported
 * index and confidence interval is reproducible from (seed, r, p, n).
 * ================================================================== */

typedef struct { uint64_t s; } Rng;

static void rng_seed(Rng* r, uint64_t seed)
{
    /* Avoid the zero fixed point of xorshift. */
    r->s = (seed == 0) ? UINT64_C(0x9E3779B97F4A7C15) : seed;
}

static uint64_t rng_next(Rng* r)
{
    uint64_t x = r->s;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    r->s = x;
    return x * UINT64_C(0x2545F4914F6CDD1D);
}

/* Uniform in [0,1). */
static double rng_uniform(Rng* r)
{
    return (double)(rng_next(r) >> 11) * (1.0 / 9007199254740992.0);
}

/* Uniform integer in [0, n). */
static size_t rng_below(Rng* r, size_t n)
{
    return (n == 0) ? 0 : (size_t)(rng_uniform(r) * (double)n) % n;
}

/* ------------------------------------------------------------------ Latin hypercube */

/*
 * Latin-hypercube sample: n points in the d-dimensional unit cube, stratified
 * so that each dimension is divided into n intervals of equal probability and
 * each interval contains exactly one point.  Column order is randomised
 * independently per dimension by a Fisher-Yates shuffle, which is what makes
 * two successive calls on the same generator independent of one another.
 *
 * Used as the base design for the Sobol' estimator when 2d exceeds the
 * dimension limit of sensitivity_generate_unit_samples(); see the comment in
 * run_sobol_parallel().
 *
 * Ref: McKay, Beckman & Conover (1979) Technometrics 21:239-245.
 */
static int latin_hypercube(size_t n, size_t d, Rng* rng, real_t* out)
{
    size_t* perm = (size_t*)malloc(n * sizeof(size_t));
    if (!perm)
        return -1;

    for (size_t k = 0; k < d; ++k) {
        for (size_t i = 0; i < n; ++i)
            perm[i] = i;

        for (size_t i = n; i-- > 1;) {
            size_t j = rng_below(rng, i + 1);
            size_t t = perm[i];
            perm[i]  = perm[j];
            perm[j]  = t;
        }

        for (size_t i = 0; i < n; ++i)
            out[i * d + k] = (real_t)(((double)perm[i] + rng_uniform(rng)) / (double)n);
    }

    free(perm);
    return 0;
}

/* Human-readable description of the design actually used, for the run banner
 * and the CSV provenance column. */
static const char* g_sobol_design = "unset";

/* ------------------------------------------------------------------ percentiles */

static int cmp_double(const void* a, const void* b)
{
    double x = *(const double*)a, y = *(const double*)b;
    return (x < y) ? -1 : (x > y) ? 1 : 0;
}

/*
 * Percentile of an *unsorted* buffer (sorted in place), linear interpolation
 * between order statistics (the "type 7" definition used by R and NumPy).
 */
static double percentile_inplace(double* v, size_t n, double q)
{
    if (n == 0) return 0.0;
    if (n == 1) return v[0];
    qsort(v, n, sizeof(double), cmp_double);
    double pos = q * (double)(n - 1);
    size_t lo  = (size_t)pos;
    size_t hi  = (lo + 1 < n) ? lo + 1 : n - 1;
    double frac = pos - (double)lo;
    return v[lo] * (1.0 - frac) + v[hi] * frac;
}

/* ==================================================================
 * Morris screening with explicit elementary effects
 *
 * Implemented locally (rather than via sensitivity_morris_analyze) so
 * that the raw elementary effects are retained: they are required for
 * the bootstrap confidence intervals on mu* requested by reviewers.
 *
 * Design: Morris (1991) one-at-a-time trajectories, mu* after
 * Campolongo, Cariboni & Saltelli (2007), doi:10.1016/j.envsoft.2006.10.004.
 *   levels p (even), step  Delta = p / (2 (p - 1))
 *   base point components drawn from {0, 1/(p-1), ..., 1 - Delta}
 *   per-parameter direction d_i in {-1, +1}, random coordinate order
 *   r trajectories  =>  r (d + 1) model evaluations
 *   EE_i = [y(z_j) - y(z_{j-1})] / (d_i Delta)      (unit-cube space)
 * ================================================================== */

typedef struct {
    real_t mu;
    real_t mu_star;
    real_t sigma;
    real_t mu_star_lo;   /* bootstrap percentile CI on mu* */
    real_t mu_star_hi;
    size_t effects_count;
} MorrisResultCI;

typedef struct {
    const real_t*             unit_row;
    real_t*                   out;
    const SensitivityProblem* problem;
    size_t                    d;
} EvalJobCtx;

static int eval_job_run(size_t job_id, void* job_ctx, void* result_ctx)
{
    (void)job_id;
    (void)result_ctx;
    EvalJobCtx* j = (EvalJobCtx*)job_ctx;

    real_t physical[MAX_PARAMS];
    if (sensitivity_scale_unit_sample(
            j->problem->parameters, j->d, j->unit_row, physical) != 0)
        return -1;

    return plant_evaluate(physical, j->d, j->out, j->problem->ctx);
}

/* Evaluate an array of unit-cube rows in parallel. */
static int evaluate_rows_parallel(
    const SensitivityProblem* problem,
    const real_t*             unit_rows,
    size_t                    row_count,
    size_t                    d,
    size_t                    workers,
    real_t*                   out_values)
{
    EvalJobCtx* ctxs = (EvalJobCtx*)calloc(row_count, sizeof(EvalJobCtx));
    BatchJob*   jobs = (BatchJob*)calloc(row_count, sizeof(BatchJob));
    int         rc   = -1;

    if (!ctxs || !jobs) goto done;

    for (size_t i = 0; i < row_count; ++i) {
        ctxs[i] = (EvalJobCtx){
            .unit_row = unit_rows + i * d,
            .out      = &out_values[i],
            .problem  = problem,
            .d        = d
        };
        jobs[i] = (BatchJob){ .run = eval_job_run, .job_ctx = &ctxs[i] };
    }

    BatchRunnerConfig cfg = { .max_workers = workers };
    rc = batch_runner_run(&cfg, jobs, row_count);

done:
    free(ctxs);
    free(jobs);
    return rc;
}

/*
 * Runs r Morris trajectories, stores every elementary effect, then
 * reduces to mu / mu* / sigma plus a bootstrap percentile interval on mu*.
 * out_effects (optional, r x d) receives the raw effects for diagnostics.
 */
static int morris_analyze_ci(
    const SensitivityProblem* problem,
    size_t                    trajectories,
    size_t                    levels,
    unsigned int              seed,
    size_t                    workers,
    size_t                    boot_reps,
    MorrisResultCI*           out_results,
    real_t*                   out_effects)
{
    size_t d = problem->parameter_count;
    size_t r = trajectories;
    size_t p = (levels < 2) ? 2 : levels;

    /* Even number of levels keeps the design balanced (Saltelli et al. 2008). */
    if (p % 2 != 0) p += 1;

    double delta = (double)p / (2.0 * (double)(p - 1));
    size_t n_base_levels = p / 2;                 /* {0, 1/(p-1), ..., 1-Delta} */
    double grid_step     = 1.0 / (double)(p - 1);

    size_t rows_per_traj = d + 1;
    size_t total_rows    = r * rows_per_traj;

    real_t* rows    = (real_t*)calloc(total_rows * d, sizeof(real_t));
    real_t* values  = (real_t*)calloc(total_rows,     sizeof(real_t));
    real_t* effects = out_effects ? out_effects
                                  : (real_t*)calloc(r * d, sizeof(real_t));
    int*    dirs    = (int*)calloc(r * d,   sizeof(int));
    size_t* orders  = (size_t*)calloc(r * d, sizeof(size_t));
    double* boot    = (double*)calloc(boot_reps ? boot_reps : 1, sizeof(double));
    int     status  = -1;

    if (!rows || !values || !effects || !dirs || !orders || !boot) goto done;

    Rng rng;
    rng_seed(&rng, (uint64_t)seed * UINT64_C(6364136223846793005) + UINT64_C(1442695040888963407));

    /* ---- Build all trajectories up front so they can be evaluated in one batch ---- */
    for (size_t t = 0; t < r; ++t) {
        real_t* z0 = rows + t * rows_per_traj * d;

        /* Base point on the coarse grid, then shifted per direction so that
         * both the point and the point + d_i*Delta remain inside [0,1]. */
        for (size_t k = 0; k < d; ++k) {
            double base = (double)rng_below(&rng, n_base_levels) * grid_step;
            int    dir  = (rng_uniform(&rng) < 0.5) ? -1 : +1;
            dirs[t * d + k] = dir;
            z0[k] = (real_t)((dir > 0) ? base : base + delta);
        }

        /* Random coordinate order (Fisher-Yates). */
        size_t* ord = orders + t * d;
        for (size_t k = 0; k < d; ++k) ord[k] = k;
        for (size_t k = d; k > 1; --k) {
            size_t j = rng_below(&rng, k);
            size_t tmp = ord[k - 1]; ord[k - 1] = ord[j]; ord[j] = tmp;
        }

        /* One-at-a-time steps. */
        for (size_t step = 0; step < d; ++step) {
            real_t* prev = rows + (t * rows_per_traj + step)     * d;
            real_t* cur  = rows + (t * rows_per_traj + step + 1) * d;
            memcpy(cur, prev, d * sizeof(real_t));
            size_t k = ord[step];
            cur[k] = (real_t)((double)cur[k] + (double)dirs[t * d + k] * delta);
            /* Guard against round-off pushing the point outside the unit cube. */
            if ((double)cur[k] < 0.0) cur[k] = REAL(0.0);
            if ((double)cur[k] > 1.0) cur[k] = REAL(1.0);
        }
    }

    if (evaluate_rows_parallel(problem, rows, total_rows, d, workers, values) != 0)
        goto done;

    /* ---- Elementary effects ---- */
    for (size_t t = 0; t < r; ++t) {
        size_t* ord = orders + t * d;
        for (size_t step = 0; step < d; ++step) {
            size_t k    = ord[step];
            double y0   = (double)values[t * rows_per_traj + step];
            double y1   = (double)values[t * rows_per_traj + step + 1];
            double sgn  = (double)dirs[t * d + k];
            effects[t * d + k] = (real_t)((y1 - y0) / (sgn * delta));
        }
    }

    /* ---- Reduce ---- */
    for (size_t k = 0; k < d; ++k) {
        double sum = 0.0, sum_abs = 0.0;
        for (size_t t = 0; t < r; ++t) {
            double e = (double)effects[t * d + k];
            sum     += e;
            sum_abs += fabs(e);
        }
        double mu      = sum / (double)r;
        double mu_star = sum_abs / (double)r;

        double var = 0.0;
        for (size_t t = 0; t < r; ++t) {
            double dv = (double)effects[t * d + k] - mu;
            var += dv * dv;
        }
        var = (r > 1) ? var / (double)(r - 1) : 0.0;

        out_results[k].mu            = (real_t)mu;
        out_results[k].mu_star       = (real_t)mu_star;
        out_results[k].sigma         = (real_t)sqrt(var);
        out_results[k].effects_count = r;

        /* Bootstrap over trajectories (the independent replicate unit). */
        if (boot_reps > 0 && r > 1) {
            Rng brng;
            rng_seed(&brng, (uint64_t)seed + UINT64_C(0x5EED0000) + (uint64_t)k);
            for (size_t b = 0; b < boot_reps; ++b) {
                double acc = 0.0;
                for (size_t t = 0; t < r; ++t) {
                    size_t idx = rng_below(&brng, r);
                    acc += fabs((double)effects[idx * d + k]);
                }
                boot[b] = acc / (double)r;
            }
            out_results[k].mu_star_lo = (real_t)percentile_inplace(boot, boot_reps, 0.025);
            out_results[k].mu_star_hi = (real_t)percentile_inplace(boot, boot_reps, 0.975);
        } else {
            out_results[k].mu_star_lo = (real_t)mu_star;
            out_results[k].mu_star_hi = (real_t)mu_star;
        }
    }

    status = 0;

done:
    free(rows);
    free(values);
    if (!out_effects) free(effects);
    free(dirs);
    free(orders);
    free(boot);
    return status;
}

/* ------------------------------------------------------------------ parallel Sobol (Saltelli 2010) */

/*
 * Sobol indices with bootstrap confidence intervals.
 *   first_order / total_order  reported (non-negative, ST clamped to 1)
 *   *_lo / *_hi                2.5 / 97.5 bootstrap percentiles
 *   total_order_raw            ST before clamping — a value > 1 is the
 *                              estimator's own warning that n is too small
 *   interaction                ST - S1, the share of variance carried by
 *                              interactions rather than the main effect
 */
typedef struct {
    real_t first_order;
    real_t total_order;
    real_t first_order_lo, first_order_hi;
    real_t total_order_lo, total_order_hi;
    real_t total_order_raw;
    real_t interaction;
} SobolResultCI;

typedef struct {
    const real_t*             unit_row;  /* pointer into unit_samples, length d */
    real_t*                   out;       /* pointer into out_values */
    const SensitivityProblem* problem;
    size_t                    d;
} SobolJobCtx;

static int sobol_job_run(size_t job_id, void* job_ctx, void* result_ctx)
{
    (void)job_id;
    (void)result_ctx;
    SobolJobCtx* j = (SobolJobCtx*)job_ctx;

    real_t physical[MAX_PARAMS];
    if (sensitivity_scale_unit_sample(
            j->problem->parameters, j->d, j->unit_row, physical) != 0)
        return -1;

    return plant_evaluate(physical, j->d, j->out, j->problem->ctx);
}

/*
 * Implements variance-based Sobol SA (Saltelli et al. 2010, doi:10.1016/j.cpc.2009.09.018).
 * Builds A, B, and AB_j matrices from a 2n quasi-random Sobol sequence, evaluates each
 * row in parallel after scaling unit [0,1] samples to physical parameter values.
 * Total evaluations: n × (d + 2).
 */
static int run_sobol_parallel(
    const SensitivityProblem* problem,
    size_t                    n,
    unsigned int              seed,
    SensitivitySamplerKind    sampler,
    size_t                    workers,
    size_t                    boot_reps,
    SobolResultCI*            out_results)
{
    size_t d          = problem->parameter_count;
    size_t total_rows = n * (d + 2);

    real_t*      unit_samples = (real_t*)calloc(total_rows * d, sizeof(real_t));
    real_t*      out_values   = (real_t*)calloc(total_rows,     sizeof(real_t));
    real_t*      sobol_2d     = (real_t*)calloc(n * 2 * d,      sizeof(real_t));
    SobolJobCtx* job_ctxs     = (SobolJobCtx*)calloc(total_rows, sizeof(SobolJobCtx));
    BatchJob*    jobs         = (BatchJob*)calloc(total_rows,     sizeof(BatchJob));
    int          status       = -1;

    if (!unit_samples || !out_values || !sobol_2d || !job_ctxs || !jobs)
        goto done;

    /* Base sample: matrices A and B, each n x d.
     *
     * The Saltelli estimator requires A and B to be INDEPENDENT samples of
     * the input distribution.  Drawing a single sequence of length 2n and
     * splitting it into halves — as an earlier version of this file did —
     * does not satisfy that requirement for a quasi-random sequence.  The two
     * halves of a Sobol sequence are strongly correlated: for the d = 16
     * problem defined in config/sensitivity_params.csv the per-column
     * correlation corr(A_k, B_k) reaches 1.0000, so replacing a column of A
     * by the corresponding column of B changes almost nothing and the
     * estimator collapses towards zero.
     *
     * Verified against the Ishigami function, whose indices are known in
     * closed form (S1 = [0.3139, 0.4424, 0], ST = [0.5576, 0.4424, 0.2437]).
     * Mean absolute error over 8 seeds, n = 256 -> 16384:
     *
     *   halves of one 2n sequence   0.31 -> 0.30   (does not converge;
     *                                               S1 for x1 estimated 0.000)
     *   column split of n x 2d seq  0.034 -> 0.0007
     *   two Sobol calls, new seed   0.042 -> 0.021 (bias floor: reseeding
     *                                               this sampler does not
     *                                               yield an independent set)
     *   two Latin-hypercube sets    0.055 -> 0.006 (converges as 1/sqrt(n))
     *
     * The column split of a single 2d-dimensional low-discrepancy sequence is
     * therefore preferred, but sensitivity_generate_unit_samples() supports
     * at most 16 dimensions, and 2d = 32 here.  The design below uses the
     * column split whenever 2d fits, and otherwise falls back to two
     * independent Latin-hypercube samples, which are stratified per dimension
     * and provably independent by construction.
     *
     * Ref: Saltelli et al. (2010) Comput. Phys. Commun. 181:259-270, sec. 2;
     *      Sobol' (2001) Math. Comput. Simul. 55:271-280;
     *      McKay, Beckman & Conover (1979) Technometrics 21:239-245.
     */
    real_t* a_block = unit_samples;         /* rows 0 .. n-1     */
    real_t* b_block = unit_samples + n * d; /* rows n .. 2n-1    */

    if (sensitivity_generate_unit_samples(sampler, n, 2 * d, seed, sobol_2d) == 0) {
        for (size_t i = 0; i < n; ++i)
            for (size_t k = 0; k < d; ++k) {
                a_block[i * d + k] = sobol_2d[i * 2 * d + k];
                b_block[i * d + k] = sobol_2d[i * 2 * d + d + k];
            }
        g_sobol_design = "column split of one n x 2d low-discrepancy sequence";
    } else {
        Rng rng;
        rng_seed(&rng, (uint64_t)seed * UINT64_C(6364136223846793005)
                     + UINT64_C(1442695040888963407));
        if (latin_hypercube(n, d, &rng, a_block) != 0) goto done;
        if (latin_hypercube(n, d, &rng, b_block) != 0) goto done;
        g_sobol_design = "two independent Latin-hypercube samples "
                         "(2d exceeds sampler dimension limit)";
    }

    /* AB_j matrices: A with column j replaced by B's column j. */
    for (size_t j = 0; j < d; ++j) {
        real_t* abj = unit_samples + (2 + j) * n * d;
        for (size_t i = 0; i < n; ++i)
            for (size_t k = 0; k < d; ++k)
                abj[i * d + k] = (k == j) ? b_block[i * d + k] : a_block[i * d + k];
    }

    /* Build parallel jobs: each scales one unit row to physical values and evaluates. */
    for (size_t i = 0; i < total_rows; ++i) {
        job_ctxs[i] = (SobolJobCtx){
            .unit_row = unit_samples + i * d,
            .out      = &out_values[i],
            .problem  = problem,
            .d        = d
        };
        jobs[i] = (BatchJob){ .run = sobol_job_run, .job_ctx = &job_ctxs[i] };
    }

    BatchRunnerConfig batch_cfg = { .max_workers = workers };
    if (batch_runner_run(&batch_cfg, jobs, total_rows) != 0)
        goto done;

    /* Saltelli estimator: mean and variance from A and B columns only. */
    real_t* ya = out_values;
    real_t* yb = out_values + n;

    real_t mean = REAL(0.0);
    for (size_t i = 0; i < n; ++i)
        mean += ya[i] + yb[i];
    mean /= (real_t)(2 * n);

    real_t variance = REAL(0.0);
    for (size_t i = 0; i < n; ++i) {
        real_t da = ya[i] - mean, db = yb[i] - mean;
        variance += da * da + db * db;
    }
    variance /= (real_t)(2 * n - 1);

    if (variance <= REAL_EPSILON) {
        fprintf(stderr,
            "Sobol: output variance ≈ 0 (mean=%.4e) — increase --samples or --days.\n",
            (double)mean);
        goto done;
    }

    /* Warn when variance is low enough that the estimator may be unreliable.
     * Relative threshold: CV < 1% of the mean suggests output saturation. */
    real_t cv = (RABS(mean) > REAL_EPSILON) ? RSQRT(variance) / RABS(mean) : REAL(0.0);
    if (cv < REAL(0.01)) {
        fprintf(stderr,
            "Sobol: low output variance (CV=%.4f, mean=%.4e) — indices may exceed 1 due to estimator instability. Increase --samples or reconsider metric.\n",
            (double)cv, (double)mean);
    }

    /* First-order S1 and total-order ST (Saltelli 2010, eqs. 21 and 23).
     * ST is theoretically in [0,1] but can exceed 1 with insufficient samples
     * or near-zero output variance; clamp and count violations for diagnostics.
     *
     * Confidence intervals: bootstrap over the n base rows.  Resampling
     * indices i in {0..n-1} keeps each triple (ya[i], yb[i], yabj[i])
     * together, so the resampled estimator stays a valid Saltelli estimator
     * on a resampled design (Archer et al. 1997; Saltelli et al. 2010).
     * Mean and variance are recomputed inside each replicate — using the
     * full-sample variance would understate the interval.
     */
    int     st_violations = 0;
    size_t* bidx  = (boot_reps > 0) ? (size_t*)calloc(n, sizeof(size_t))         : NULL;
    double* b_s1  = (boot_reps > 0) ? (double*)calloc(boot_reps, sizeof(double)) : NULL;
    double* b_st  = (boot_reps > 0) ? (double*)calloc(boot_reps, sizeof(double)) : NULL;
    double* b_tmp = (boot_reps > 0) ? (double*)calloc(boot_reps, sizeof(double)) : NULL;

    for (size_t j = 0; j < d; ++j) {
        real_t* yabj = out_values + (2 + j) * n;
        real_t  s1   = REAL(0.0), st = REAL(0.0);
        for (size_t i = 0; i < n; ++i) {
            s1 += yb[i] * (yabj[i] - ya[i]);
            real_t diff = ya[i] - yabj[i];
            st += diff * diff;
        }
        real_t st_raw = st / (REAL(2.0) * (real_t)n * variance);
        if (st_raw > REAL(1.0))
            st_violations++;

        real_t s1_val = RMAX(s1 / ((real_t)n * variance), REAL(0.0));
        real_t st_val = RMIN(RMAX(st_raw, REAL(0.0)), REAL(1.0));

        out_results[j].first_order     = s1_val;
        out_results[j].total_order     = st_val;
        out_results[j].total_order_raw = st_raw;
        out_results[j].interaction     = RMAX(st_val - s1_val, REAL(0.0));

        if (boot_reps > 0 && bidx && b_s1 && b_st && b_tmp) {
            Rng brng;
            rng_seed(&brng, (uint64_t)seed + UINT64_C(0xB007) + (uint64_t)j * UINT64_C(7919));
            for (size_t b = 0; b < boot_reps; ++b) {
                for (size_t i = 0; i < n; ++i)
                    bidx[i] = rng_below(&brng, n);

                double bmean = 0.0;
                for (size_t i = 0; i < n; ++i)
                    bmean += (double)ya[bidx[i]] + (double)yb[bidx[i]];
                bmean /= (double)(2 * n);

                double bvar = 0.0;
                for (size_t i = 0; i < n; ++i) {
                    double da = (double)ya[bidx[i]] - bmean;
                    double db = (double)yb[bidx[i]] - bmean;
                    bvar += da * da + db * db;
                }
                bvar /= (double)(2 * n - 1);

                if (!(bvar > 0.0)) { b_s1[b] = (double)s1_val; b_st[b] = (double)st_val; continue; }

                double acc1 = 0.0, acct = 0.0;
                for (size_t i = 0; i < n; ++i) {
                    size_t k = bidx[i];
                    acc1 += (double)yb[k] * ((double)yabj[k] - (double)ya[k]);
                    double dv = (double)ya[k] - (double)yabj[k];
                    acct += dv * dv;
                }
                b_s1[b] = fmax(acc1 / ((double)n * bvar), 0.0);
                b_st[b] = fmin(fmax(acct / (2.0 * (double)n * bvar), 0.0), 1.0);
            }

            memcpy(b_tmp, b_s1, boot_reps * sizeof(double));
            out_results[j].first_order_lo = (real_t)percentile_inplace(b_tmp, boot_reps, 0.025);
            memcpy(b_tmp, b_s1, boot_reps * sizeof(double));
            out_results[j].first_order_hi = (real_t)percentile_inplace(b_tmp, boot_reps, 0.975);
            memcpy(b_tmp, b_st, boot_reps * sizeof(double));
            out_results[j].total_order_lo = (real_t)percentile_inplace(b_tmp, boot_reps, 0.025);
            memcpy(b_tmp, b_st, boot_reps * sizeof(double));
            out_results[j].total_order_hi = (real_t)percentile_inplace(b_tmp, boot_reps, 0.975);
        } else {
            out_results[j].first_order_lo = s1_val;
            out_results[j].first_order_hi = s1_val;
            out_results[j].total_order_lo = st_val;
            out_results[j].total_order_hi = st_val;
        }
    }

    free(bidx);
    free(b_s1);
    free(b_st);
    free(b_tmp);

    if (st_violations > 0) {
        fprintf(stderr,
            "Sobol: %d ST index(es) exceeded 1 before clamping — results unreliable, increase --samples.\n",
            st_violations);
    }

    status = 0;

done:
    free(unit_samples);
    free(out_values);
    free(sobol_2d);
    free(job_ctxs);
    free(jobs);
    return status;
}

/* ------------------------------------------------------------------ helpers */

static void print_usage(const char* argv0)
{
    printf(
        "Usage: %s [options]\n"
        "\n"
        "Options:\n"
        "  --config <path>         Base simulation config  [config/default.config]\n"
        "  --params <path>         Sensitivity params CSV  [config/sensitivity_params.csv]\n"
        "  --analysis <kind>       morris | sobol | both   [both]\n"
        "  --metric <kind>         rgr | biomass | sucrose | starch | carb  [rgr]\n"
        "  --photoperiod <h>       Fixed photoperiod hours [12]\n"
        "  --all-photoperiods      Run for 4, 6, 8, 12 h and write per-photoperiod CSVs\n"
        "  --days <N>              Simulation days         [30]\n"
        "  --samples <N>           Sobol sample count      [512]\n"
        "  --trajectories <N>      Morris trajectories     [30]\n"
        "  --levels <N>            Morris levels           [6]\n"
        "  --workers <N>           Parallel workers        [4]\n"
        "  --seed <N>              Random seed             [0]\n"
        "  --sampler sobol|halton  Quasi-random sampler    [sobol]\n"
        "  --bootstrap <n>         Bootstrap replicates    [1000, 0 disables]\n"
        "  --out-morris <path>     Morris result CSV       [sensitivity_morris.csv]\n"
        "  --out-sobol <path>      Sobol result CSV        [sensitivity_sobol.csv]\n"
        "\n"
        "Notes:\n"
        "  - Sobol max dimension: 16 (solver library limit).\n"
        "  - Sobol total evaluations: samples * (params + 2), run in parallel.\n"
        "  - Morris total evaluations: 2 * trajectories * params (sequential).\n",
        argv0);
}

static const char* metric_label(MetricKind m)
{
    switch (m) {
        case METRIC_RGR:     return "RGR [day^-1]";
        case METRIC_BIOMASS: return "Final biomass [gFW]";
        case METRIC_SUCROSE: return "Final sucrose [gFW gDW^-1]";
        case METRIC_STARCH:  return "Final starch [gFW gDW^-1]";
        case METRIC_CARB:    return "J_carb = 0.5*(sucrose/max_sucrose + starch/max_starch)";
        default:             return "?";
    }
}

static void elapsed_str(double sec, char* buf, size_t len)
{
    int h = (int)(sec / 3600.0);
    int m = (int)((sec - h * 3600.0) / 60.0);
    int s = (int)(sec - h * 3600.0 - m * 60.0);
    if (h > 0)
        snprintf(buf, len, "%dh %dm %ds", h, m, s);
    else if (m > 0)
        snprintf(buf, len, "%dm %ds", m, s);
    else
        snprintf(buf, len, "%ds", s);
}

/* ------------------------------------------------------------------ CSV writers
 *
 * Written locally rather than via sensitivity_write_*_csv so that the
 * confidence-interval columns are included.  Column order is kept
 * backward-compatible with the previous files: the new columns are appended,
 * so existing plotting scripts keep working.
 */

static void write_morris_ci_csv(
    FILE*                       f,
    const SensitivityParameter* params,
    const MorrisResultCI*       res,
    size_t                      count,
    size_t                      levels,
    size_t                      boot_reps)
{
    fprintf(f, "name,role,scale,lower,upper,mu,mu_star,sigma,"
               "mu_star_lo,mu_star_hi,effects_count,levels,bootstrap_reps\n");
    for (size_t i = 0; i < count; ++i) {
        fprintf(f, "%s,%s,%s,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%zu,%zu,%zu\n",
            params[i].name,
            sensitivity_role_string(params[i].role),
            sensitivity_scale_string(params[i].scale),
            (double)params[i].lower,
            (double)params[i].upper,
            (double)res[i].mu,
            (double)res[i].mu_star,
            (double)res[i].sigma,
            (double)res[i].mu_star_lo,
            (double)res[i].mu_star_hi,
            res[i].effects_count,
            levels,
            boot_reps);
    }
}

static void write_sobol_ci_csv(
    FILE*                       f,
    const SensitivityParameter* params,
    const SobolResultCI*        res,
    size_t                      count,
    size_t                      sample_count,
    size_t                      boot_reps)
{
    fprintf(f, "name,role,scale,lower,upper,S1,ST,sample_count,"
               "S1_lo,S1_hi,ST_lo,ST_hi,ST_raw,interaction,bootstrap_reps,base_design\n");
    for (size_t i = 0; i < count; ++i) {
        fprintf(f, "%s,%s,%s,%.17g,%.17g,%.17g,%.17g,%zu,"
                   "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%zu",
            params[i].name,
            sensitivity_role_string(params[i].role),
            sensitivity_scale_string(params[i].scale),
            (double)params[i].lower,
            (double)params[i].upper,
            (double)res[i].first_order,
            (double)res[i].total_order,
            sample_count,
            (double)res[i].first_order_lo,
            (double)res[i].first_order_hi,
            (double)res[i].total_order_lo,
            (double)res[i].total_order_hi,
            (double)res[i].total_order_raw,
            (double)res[i].interaction,
            boot_reps);
        fprintf(f, ",\"%s\"\n", g_sobol_design);
    }
}

/* Raw elementary effects, one row per (trajectory, parameter). */
static void write_morris_effects_csv(
    FILE*                       f,
    const SensitivityParameter* params,
    const real_t*               effects,
    size_t                      trajectories,
    size_t                      d)
{
    fprintf(f, "trajectory,name,elementary_effect\n");
    for (size_t t = 0; t < trajectories; ++t)
        for (size_t k = 0; k < d; ++k)
            fprintf(f, "%zu,%s,%.17g\n",
                t, params[k].name, (double)effects[t * d + k]);
}

/* Build output path: inject photoperiod suffix before extension. */
static void build_out_path(const char* base, int ph_h, char* dst, size_t dsz)
{
    if (ph_h <= 0) {
        snprintf(dst, dsz, "%s", base);
        return;
    }
    const char* ext = strrchr(base, '.');
    if (ext && strcmp(ext, ".csv") == 0) {
        size_t prefix = (size_t)(ext - base);
        snprintf(dst, dsz, "%.*s_%dh.csv", (int)prefix, base, ph_h);
    } else {
        snprintf(dst, dsz, "%s_%dh.csv", base, ph_h);
    }
}

/* ------------------------------------------------------------------ single run */

static int run_for_photoperiod(
    int                         ph_h,
    const SensitivityParameter* params,
    size_t                      param_count,
    const PlantSensCtx*         base_ctx,
    AnalysisKind                analysis,
    int                         morris_traj,
    int                         morris_levels,
    size_t                      sobol_n,
    size_t                      boot_reps,
    unsigned int                seed,
    SensitivitySamplerKind      sampler,
    size_t                      workers,
    const char*                 out_morris,
    const char*                 out_sobol)
{
    PlantSensCtx ctx = *base_ctx;
    if (ph_h > 0)
        ctx.base_input.photo.photoperiod = (real_t)ph_h;

    SensitivityProblem problem = {
        .parameters      = params,
        .parameter_count = param_count,
        .evaluate        = plant_evaluate,
        .ctx             = &ctx
    };

    const char* tag = (ph_h > 0) ? "" : "";
    char tag_buf[32] = "";
    if (ph_h > 0)
        snprintf(tag_buf, sizeof(tag_buf), " [photoperiod=%dh]", ph_h);
    tag = tag_buf;

    int any_error = 0;

    /* ---- Morris ---- */
    if (analysis & ANALYSIS_MORRIS) {
        printf("\n=== Morris screening%s ===\n", tag);
        printf("  trajectories=%d  levels=%d  evaluations=%zu\n",
            morris_traj, morris_levels,
            (size_t)(2 * morris_traj) * param_count);

        atomic_store(&ctx.n_eval, 0UL);
        atomic_store(&ctx.n_diverged, 0UL);
            atomic_store(&ctx.n_bad_endpoint, 0UL);

        MorrisResultCI* morris_res = (MorrisResultCI*)calloc(param_count, sizeof(MorrisResultCI));
        real_t*         effects    = (real_t*)calloc((size_t)morris_traj * param_count, sizeof(real_t));
        if (!morris_res || !effects) {
            fprintf(stderr, "OOM\n");
            free(morris_res); free(effects);
            return -1;
        }

        time_t t0 = time(NULL);
        int rc = morris_analyze_ci(&problem,
                    (size_t)morris_traj, (size_t)morris_levels,
                    seed, workers, boot_reps,
                    morris_res, effects);
        double elapsed = difftime(time(NULL), t0);
        char etime[32];
        elapsed_str(elapsed, etime, sizeof(etime));

        if (rc != 0) {
            fprintf(stderr, "Morris analysis failed.\n");
            free(morris_res);
            free(effects);
            any_error = 1;
        } else {
            printf("\n  %-44s  %10s  %10s  %10s  %21s\n",
                "Parameter", "mu", "mu*", "sigma", "mu* 95% CI");
            printf("  %s\n", "-----------------------------------------------------------------------------------------------------------------");

            /* Sort copy by mu* descending for readability */
            size_t* order = (size_t*)calloc(param_count, sizeof(size_t));
            if (order) {
                for (size_t i = 0; i < param_count; ++i) order[i] = i;
                /* Bubble sort (small n) */
                for (size_t i = 0; i < param_count; ++i)
                    for (size_t j = i + 1; j < param_count; ++j)
                        if ((double)morris_res[order[j]].mu_star > (double)morris_res[order[i]].mu_star) {
                            size_t tmp = order[i]; order[i] = order[j]; order[j] = tmp;
                        }

                for (size_t r = 0; r < param_count; ++r) {
                    size_t i = order[r];
                    printf("  %-44s  %10.4f  %10.4f  %10.4f  [%8.4f, %8.4f]\n",
                        params[i].name,
                        (double)morris_res[i].mu,
                        (double)morris_res[i].mu_star,
                        (double)morris_res[i].sigma,
                        (double)morris_res[i].mu_star_lo,
                        (double)morris_res[i].mu_star_hi);
                }
                free(order);
            }

            unsigned long ne = atomic_load(&ctx.n_eval);
            unsigned long nd = atomic_load(&ctx.n_diverged);
            unsigned long nb = atomic_load(&ctx.n_bad_endpoint);
            printf("  diverged evaluations: %lu of %lu (%.2f%%)"
                   "; unusable endpoint: %lu (%.2f%%)\n",
                nd, ne, ne ? 100.0 * (double)nd / (double)ne : 0.0,
                nb, ne ? 100.0 * (double)nb / (double)ne : 0.0);
            printf("\n  Completed in %s\n", etime);

            char path[512];
            build_out_path(out_morris, ph_h, path, sizeof(path));
            FILE* f = fopen(path, "w");
            if (f) {
                write_morris_ci_csv(f, params, morris_res, param_count,
                    (size_t)morris_levels, boot_reps);
                fclose(f);
                printf("  → %s\n", path);
            } else {
                fprintf(stderr, "  Could not write %s\n", path);
            }

            /* Raw elementary effects alongside the summary, for the
             * convergence/CI diagnostics reported in the appendix. */
            char epath[600];
            snprintf(epath, sizeof(epath), "%.*s_effects%s",
                (int)(strlen(path) - 4), path, ".csv");
            FILE* fe = fopen(epath, "w");
            if (fe) {
                write_morris_effects_csv(fe, params, effects,
                    (size_t)morris_traj, param_count);
                fclose(fe);
                printf("  → %s\n", epath);
            }
            free(morris_res);
            free(effects);
        }
    }

    /* ---- Sobol ---- */
    if (analysis & ANALYSIS_SOBOL) {
        if (param_count > MAX_PARAMS) {
            fprintf(stderr,
                "Sobol: %zu parameters exceeds solver dimension limit %d — skipping.\n",
                param_count, MAX_PARAMS);
        } else {
            size_t total_evals = sobol_n * (param_count + 2);
            printf("\n=== Sobol variance decomposition%s ===\n", tag);
            atomic_store(&ctx.n_eval, 0UL);
            atomic_store(&ctx.n_diverged, 0UL);
            atomic_store(&ctx.n_bad_endpoint, 0UL);
            printf("  samples=%zu  workers=%zu  total_evaluations=%zu\n",
                sobol_n, workers, total_evals);
            printf("  sampler=%s  seed=%u\n",
                (sampler == SENSITIVITY_SAMPLER_SOBOL_SEQUENCE) ? "sobol"
                : (sampler == SENSITIVITY_SAMPLER_LOW_DISCREPANCY) ? "halton" : "monte-carlo",
                seed);

            SobolResultCI* sobol_res = (SobolResultCI*)calloc(param_count, sizeof(SobolResultCI));
            if (!sobol_res) { fprintf(stderr, "OOM\n"); return -1; }

            time_t t0 = time(NULL);
            int rc = run_sobol_parallel(&problem, sobol_n, seed, sampler, workers,
                        boot_reps, sobol_res);
            double elapsed = difftime(time(NULL), t0);
            printf("  base design: %s\n", g_sobol_design);
            char etime[32];
            elapsed_str(elapsed, etime, sizeof(etime));

            if (rc != 0) {
                fprintf(stderr, "Sobol analysis failed.\n");
                any_error = 1;
            } else {
                printf("\n  %-44s  %10s  %19s  %10s  %19s  %10s\n",
                    "Parameter", "S1", "S1 95% CI", "ST", "ST 95% CI", "ST-S1");
                printf("  %s\n",
                    "--------------------------------------------------------------------------------------------------------------------------------");

                size_t* order = (size_t*)calloc(param_count, sizeof(size_t));
                if (order) {
                    for (size_t i = 0; i < param_count; ++i) order[i] = i;
                    for (size_t i = 0; i < param_count; ++i)
                        for (size_t j = i + 1; j < param_count; ++j)
                            if ((double)sobol_res[order[j]].total_order > (double)sobol_res[order[i]].total_order) {
                                size_t tmp = order[i]; order[i] = order[j]; order[j] = tmp;
                            }

                    for (size_t r = 0; r < param_count; ++r) {
                        size_t i = order[r];
                        printf("  %-44s  %10.4f  [%7.4f, %7.4f]  %10.4f  [%7.4f, %7.4f]  %10.4f\n",
                            params[i].name,
                            (double)sobol_res[i].first_order,
                            (double)sobol_res[i].first_order_lo,
                            (double)sobol_res[i].first_order_hi,
                            (double)sobol_res[i].total_order,
                            (double)sobol_res[i].total_order_lo,
                            (double)sobol_res[i].total_order_hi,
                            (double)sobol_res[i].interaction);
                    }
                    free(order);
                }

                unsigned long ne = atomic_load(&ctx.n_eval);
                unsigned long nd = atomic_load(&ctx.n_diverged);
                unsigned long nb = atomic_load(&ctx.n_bad_endpoint);
                printf("  diverged evaluations: %lu of %lu (%.2f%%)"
                       "; unusable endpoint: %lu (%.2f%%)\n",
                    nd, ne, ne ? 100.0 * (double)nd / (double)ne : 0.0,
                    nb, ne ? 100.0 * (double)nb / (double)ne : 0.0);
                printf("\n  Completed in %s\n", etime);

                char path[512];
                build_out_path(out_sobol, ph_h, path, sizeof(path));
                FILE* f = fopen(path, "w");
                if (f) {
                    write_sobol_ci_csv(f, params, sobol_res, param_count, sobol_n, boot_reps);
                    fclose(f);
                    printf("  → %s\n", path);
                } else {
                    fprintf(stderr, "  Could not write %s\n", path);
                }
            }
            free(sobol_res);
        }
    }

    return any_error ? -1 : 0;
}

/* ------------------------------------------------------------------ main */

int main(int argc, char** argv)
{
    const char* config_path  = "config/default.config";
    const char* params_path  = "config/sensitivity_params.csv";
    const char* out_morris   = "results/sensitivity/morris/morris.csv";
    const char* out_sobol    = "results/sensitivity/sobol/sobol.csv";
    const char* analysis_str = "both";
    const char* metric_str   = "rgr";
    const char* sampler_str  = "sobol";

    int    ph_h             = 12;
    int    all_photoperiods = 0;
    int    days_override    = 0;
    size_t sobol_n          = 512;
    size_t boot_reps        = 1000;
    int    morris_traj      = 30;
    int    morris_levels    = 6;
    size_t workers          = 4;
    unsigned int seed       = 0;

    for (int i = 1; i < argc; ++i) {
#define NEXTARG(dst)                                                   \
    do {                                                               \
        if (i + 1 >= argc) {                                           \
            fprintf(stderr, "Missing value for %s\n", argv[i]);       \
            return 1;                                                  \
        }                                                              \
        (dst) = argv[++i];                                             \
    } while (0)
#define NEXTINT(dst)                                                   \
    do {                                                               \
        const char* _s; NEXTARG(_s);                                   \
        (dst) = (int)strtol(_s, NULL, 10);                            \
    } while (0)
#define NEXTSZ(dst)                                                    \
    do {                                                               \
        const char* _s; NEXTARG(_s);                                   \
        (dst) = (size_t)strtoul(_s, NULL, 10);                        \
    } while (0)

        if      (strcmp(argv[i], "--config")           == 0) NEXTARG(config_path);
        else if (strcmp(argv[i], "--params")           == 0) NEXTARG(params_path);
        else if (strcmp(argv[i], "--out-morris")       == 0) NEXTARG(out_morris);
        else if (strcmp(argv[i], "--out-sobol")        == 0) NEXTARG(out_sobol);
        else if (strcmp(argv[i], "--analysis")         == 0) NEXTARG(analysis_str);
        else if (strcmp(argv[i], "--metric")           == 0) NEXTARG(metric_str);
        else if (strcmp(argv[i], "--sampler")          == 0) NEXTARG(sampler_str);
        else if (strcmp(argv[i], "--photoperiod")      == 0) NEXTINT(ph_h);
        else if (strcmp(argv[i], "--days")             == 0) NEXTINT(days_override);
        else if (strcmp(argv[i], "--trajectories")     == 0) NEXTINT(morris_traj);
        else if (strcmp(argv[i], "--levels")           == 0) NEXTINT(morris_levels);
        else if (strcmp(argv[i], "--samples")          == 0) NEXTSZ(sobol_n);
        else if (strcmp(argv[i], "--bootstrap")        == 0) NEXTSZ(boot_reps);
        else if (strcmp(argv[i], "--workers")          == 0) NEXTSZ(workers);
        else if (strcmp(argv[i], "--seed")             == 0) { const char* _s; NEXTARG(_s); seed = (unsigned int)strtoul(_s, NULL, 10); }
        else if (strcmp(argv[i], "--all-photoperiods") == 0) all_photoperiods = 1;
        else if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
            print_usage(argv[0]);
            return 0;
        } else {
            fprintf(stderr, "Unknown argument: %s\n", argv[i]);
            print_usage(argv[0]);
            return 1;
        }
#undef NEXTARG
#undef NEXTINT
#undef NEXTSZ
    }

    /* Parse analysis kind */
    AnalysisKind analysis;
    if      (strcmp(analysis_str, "morris") == 0) analysis = ANALYSIS_MORRIS;
    else if (strcmp(analysis_str, "sobol")  == 0) analysis = ANALYSIS_SOBOL;
    else if (strcmp(analysis_str, "both")   == 0) analysis = ANALYSIS_BOTH;
    else {
        fprintf(stderr, "Unknown analysis kind: %s (use morris|sobol|both)\n", analysis_str);
        return 1;
    }

    /* Parse metric */
    MetricKind metric;
    if      (strcmp(metric_str, "rgr")     == 0) metric = METRIC_RGR;
    else if (strcmp(metric_str, "biomass") == 0) metric = METRIC_BIOMASS;
    else if (strcmp(metric_str, "sucrose") == 0) metric = METRIC_SUCROSE;
    else if (strcmp(metric_str, "starch")  == 0) metric = METRIC_STARCH;
    else if (strcmp(metric_str, "carb")    == 0) metric = METRIC_CARB;
    else {
        fprintf(stderr, "Unknown metric: %s (use rgr|biomass|sucrose|starch|carb)\n", metric_str);
        return 1;
    }

    /* Parse sampler */
    SensitivitySamplerKind sampler;
    if      (strcmp(sampler_str, "sobol")  == 0) sampler = SENSITIVITY_SAMPLER_SOBOL_SEQUENCE;
    else if (strcmp(sampler_str, "halton") == 0) sampler = SENSITIVITY_SAMPLER_LOW_DISCREPANCY;
    else if (strcmp(sampler_str, "mc")     == 0) sampler = SENSITIVITY_SAMPLER_MONTE_CARLO;
    else {
        fprintf(stderr, "Unknown sampler: %s (use sobol|halton|mc)\n", sampler_str);
        return 1;
    }

    /* Clamp workers */
    if (workers < 1)           workers = 1;
    if (workers > MAX_WORKERS) workers = MAX_WORKERS;

    /* ---- Load base config ---- */
    SimulationConfig cfg   = simulation_default_config();
    Input            input = generate_input();

    int force_rk4 = 0;
    for (int i = 1; i < argc; ++i)
        if (strcmp(argv[i], "--solver") == 0 && i + 1 < argc
            && strcmp(argv[i + 1], "rk4") == 0) force_rk4 = 1;

    if (config_load_file(config_path, &cfg, &input) != 0) {
        fprintf(stderr, "Failed to load config: %s\n", config_path);
        return 1;
    }

    cfg.gui_enabled = 0;
    cfg.write_csv   = 0;
    /* Solver choice, selectable because it decides what the indices measure.
     * The previous unconditional SOLVER_RK4 was justified as preventing
     * step-size collapse at the design extremes. Measured, it does the opposite:
     * at the production dt = 0.1 h the fixed step is unstable over much of the
     * design, the instability lands in the metric as a clamped RGR of -1.094
     * against a physical range of 0.07-0.45, and the resulting indicator
     * variable carries a large share of the decomposed variance. Every failing
     * point tested recovers at dt/2 and under adaptive ODE45, converging on the
     * same RGR, so the failures were a solver artefact and not a domain exit.
     * Default is now the adaptive solver from the configuration; --solver rk4
     * restores the old behaviour for comparison. */
    if (force_rk4)
        cfg.solver_type = SOLVER_RK4;
    if (days_override > 0)
        cfg.days = days_override;

    /* ---- Load sensitivity parameters ---- */
    SensitivityParameter params[MAX_PARAMS];
    size_t               param_count = 0;

    if (sensitivity_load_parameter_csv(params_path, params, MAX_PARAMS, &param_count) != 0) {
        fprintf(stderr, "Failed to load sensitivity params: %s (capacity=%d)\n",
            params_path, MAX_PARAMS);
        return 1;
    }
    if (param_count == 0) {
        fprintf(stderr, "No enabled parameters found in %s\n", params_path);
        return 1;
    }

    /* ---- Resolve Input field offsets ---- */
    size_t offsets[MAX_PARAMS];
    for (size_t i = 0; i < param_count; ++i) {
        if (config_get_input_field_offset(params[i].name, &offsets[i]) != 0) {
            fprintf(stderr, "Unknown input field: '%s'\n", params[i].name);
            return 1;
        }
    }

    /* ---- Print run summary ---- */
    printf("=== Plant Model Sensitivity Analysis ===\n");
    printf("  Config:    %s\n", config_path);
    printf("  Params:    %s  (%zu active)\n", params_path, param_count);
    printf("  Metric:    %s\n", metric_label(metric));
    printf("  Days:      %d\n", cfg.days);
    printf("  Analysis:  %s\n", analysis_str);
    if (analysis & ANALYSIS_MORRIS)
        printf("  Morris:    %d trajectories × %d levels\n", morris_traj, morris_levels);
    if (analysis & ANALYSIS_SOBOL)
        printf("  Sobol:     %zu samples × %zu workers  (total evals: %zu)\n",
            sobol_n, workers, sobol_n * (param_count + 2));
    printf("  Bootstrap: %zu replicates  (95%% percentile intervals)\n", boot_reps);
    printf("  Seed:      %u\n", seed);

    printf("\n  Parameters:\n");
    for (size_t i = 0; i < param_count; ++i) {
        printf("    [%2zu] %-44s  [%.4g, %.4g]  %s  %s\n",
            i + 1,
            params[i].name,
            (double)params[i].lower,
            (double)params[i].upper,
            sensitivity_scale_string(params[i].scale),
            sensitivity_role_string(params[i].role));
    }

    /* ---- Build context ---- */
    PlantSensCtx base_ctx;
    base_ctx.base_input = input;
    base_ctx.cfg        = cfg;
    base_ctx.metric     = metric;
    for (size_t i = 0; i < param_count; ++i)
        base_ctx.offsets[i] = offsets[i];

    /* ---- Run analysis ---- */
    static const Scenario PHOTOPERIODS[] = {
        { 4,  0.0680 },
        { 6,  0.1135 },
        { 8,  0.1708 },
        { 12, 0.2600 }
    };

    int rc = 0;

    if (all_photoperiods) {
        size_t n_ph = sizeof(PHOTOPERIODS) / sizeof(PHOTOPERIODS[0]);
        for (size_t p = 0; p < n_ph; ++p) {
            rc |= run_for_photoperiod(
                PHOTOPERIODS[p].photoperiod_h,
                params, param_count,
                &base_ctx,
                analysis,
                morris_traj, morris_levels,
                sobol_n, boot_reps, seed, sampler, workers,
                out_morris, out_sobol);
        }
    } else {
        rc = run_for_photoperiod(
            ph_h,
            params, param_count,
            &base_ctx,
            analysis,
            morris_traj, morris_levels,
            sobol_n, boot_reps, seed, sampler, workers,
            out_morris, out_sobol);
    }

    printf("\n=== Done ===\n");
    return rc == 0 ? 0 : 1;
}

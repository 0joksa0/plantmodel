#define SOLVER_REAL_AS double
#include "config/config.h"
#include "model/input.h"
#include "model/model.h"
#include "simulation/simulation.h"

#include <framework/batch_runner.h>
#include <framework/sensitivity.h>

#include <math.h>
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

/* ------------------------------------------------------------------ parallel Sobol (Saltelli 2010) */

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
    SobolResult*              out_results)
{
    size_t d          = problem->parameter_count;
    size_t total_rows = n * (d + 2);

    real_t*      unit_samples = (real_t*)calloc(total_rows * d, sizeof(real_t));
    real_t*      out_values   = (real_t*)calloc(total_rows,     sizeof(real_t));
    real_t*      sobol_2n     = (real_t*)calloc(2 * n * d,      sizeof(real_t));
    SobolJobCtx* job_ctxs     = (SobolJobCtx*)calloc(total_rows, sizeof(SobolJobCtx));
    BatchJob*    jobs         = (BatchJob*)calloc(total_rows,     sizeof(BatchJob));
    int          status       = -1;

    if (!unit_samples || !out_values || !sobol_2n || !job_ctxs || !jobs)
        goto done;

    /* Generate 2n×d quasi-random unit points; split into A (first half) and B (second half). */
    if (sensitivity_generate_unit_samples(sampler, 2 * n, d, seed, sobol_2n) != 0)
        goto done;

    real_t* a_block = unit_samples;         /* rows 0 .. n-1 */
    real_t* b_block = unit_samples + n * d; /* rows n .. 2n-1 */

    for (size_t i = 0; i < n; ++i)
        for (size_t k = 0; k < d; ++k) {
            a_block[i * d + k] = sobol_2n[i * d + k];
            b_block[i * d + k] = sobol_2n[(n + i) * d + k];
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
     * or near-zero output variance; clamp and count violations for diagnostics. */
    int st_violations = 0;
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
        out_results[j].first_order = RMAX(s1 / ((real_t)n * variance), REAL(0.0));
        out_results[j].total_order = RMIN(RMAX(st_raw, REAL(0.0)), REAL(1.0));
    }
    if (st_violations > 0) {
        fprintf(stderr,
            "Sobol: %d ST index(es) exceeded 1 before clamping — results unreliable, increase --samples.\n",
            st_violations);
    }

    status = 0;

done:
    free(unit_samples);
    free(out_values);
    free(sobol_2n);
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

        MorrisConfig morris_cfg = {
            .levels       = (size_t)morris_levels,
            .trajectories = (size_t)morris_traj,
            .seed         = seed
        };
        MorrisResult* morris_res = (MorrisResult*)calloc(param_count, sizeof(MorrisResult));
        if (!morris_res) { fprintf(stderr, "OOM\n"); return -1; }

        time_t t0 = time(NULL);
        int rc = sensitivity_morris_analyze(&problem, &morris_cfg, morris_res);
        double elapsed = difftime(time(NULL), t0);
        char etime[32];
        elapsed_str(elapsed, etime, sizeof(etime));

        if (rc != 0) {
            fprintf(stderr, "Morris analysis failed.\n");
            free(morris_res);
            any_error = 1;
        } else {
            printf("\n  %-44s  %10s  %10s  %10s\n", "Parameter", "mu", "mu*", "sigma");
            printf("  %s\n", "---------------------------------------------------------------------------------------------");

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
                    printf("  %-44s  %10.4f  %10.4f  %10.4f\n",
                        params[i].name,
                        (double)morris_res[i].mu,
                        (double)morris_res[i].mu_star,
                        (double)morris_res[i].sigma);
                }
                free(order);
            }

            printf("\n  Completed in %s\n", etime);

            char path[512];
            build_out_path(out_morris, ph_h, path, sizeof(path));
            FILE* f = fopen(path, "w");
            if (f) {
                sensitivity_write_morris_csv(f, params, morris_res, param_count);
                fclose(f);
                printf("  → %s\n", path);
            } else {
                fprintf(stderr, "  Could not write %s\n", path);
            }
        }
        free(morris_res);
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
            printf("  samples=%zu  workers=%zu  total_evaluations=%zu\n",
                sobol_n, workers, total_evals);
            printf("  sampler=%s  seed=%u\n",
                (sampler == SENSITIVITY_SAMPLER_SOBOL_SEQUENCE) ? "sobol"
                : (sampler == SENSITIVITY_SAMPLER_LOW_DISCREPANCY) ? "halton" : "monte-carlo",
                seed);

            SobolResult* sobol_res = (SobolResult*)calloc(param_count, sizeof(SobolResult));
            if (!sobol_res) { fprintf(stderr, "OOM\n"); return -1; }

            time_t t0 = time(NULL);
            int rc = run_sobol_parallel(&problem, sobol_n, seed, sampler, workers, sobol_res);
            double elapsed = difftime(time(NULL), t0);
            char etime[32];
            elapsed_str(elapsed, etime, sizeof(etime));

            if (rc != 0) {
                fprintf(stderr, "Sobol analysis failed.\n");
                any_error = 1;
            } else {
                printf("\n  %-44s  %10s  %10s\n", "Parameter", "S1", "ST");
                printf("  %s\n", "--------------------------------------------------------------------");

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
                        printf("  %-44s  %10.4f  %10.4f\n",
                            params[i].name,
                            (double)sobol_res[i].first_order,
                            (double)sobol_res[i].total_order);
                    }
                    free(order);
                }

                printf("\n  Completed in %s\n", etime);

                char path[512];
                build_out_path(out_sobol, ph_h, path, sizeof(path));
                FILE* f = fopen(path, "w");
                if (f) {
                    sensitivity_write_sobol_csv(f, params, sobol_res, param_count, sobol_n);
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

    if (config_load_file(config_path, &cfg, &input) != 0) {
        fprintf(stderr, "Failed to load config: %s\n", config_path);
        return 1;
    }

    cfg.gui_enabled = 0;
    cfg.write_csv   = 0;
    cfg.solver_type = SOLVER_RK4;   /* fixed step: prevents step-size collapse at Sobol parameter extremes */
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
                sobol_n, seed, sampler, workers,
                out_morris, out_sobol);
        }
    } else {
        rc = run_for_photoperiod(
            ph_h,
            params, param_count,
            &base_ctx,
            analysis,
            morris_traj, morris_levels,
            sobol_n, seed, sampler, workers,
            out_morris, out_sobol);
    }

    printf("\n=== Done ===\n");
    return rc == 0 ? 0 : 1;
}

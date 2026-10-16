/*
 * horizon_study.c — sensitivity of the reported RGR to the simulation horizon
 *                   and to excluding an initial spin-up
 *
 * Motivation
 * ----------
 * Four places in this repository disagree about how long a validation run is
 * and how the growth rate is extracted from it:
 *
 *   config/default.config       simulation.days = 35
 *   tools/validation_runner.c   days = 30, divisor (days - 1) = 29
 *   tests/test.c                config.days = 30, divisor 29, different
 *                               lambda_sb values than validation_runner
 *   tools/carbon_balance.c      days = 35
 *
 * This tool measures what that choice is worth in RGR, and what lambda_sb
 * would have to be under each candidate convention to reproduce the published
 * Sulpice et al. (2014) growth rates.
 *
 * The divisor question
 * --------------------
 * compute_rgr() (src/simulation/simulation.c:487) is
 *
 *     compute_rgr(B0, B1, D) = (log B1 - log B0) / D
 *
 * i.e. a plain log ratio over a caller-supplied duration.  The callers in
 * validation_runner.c, ablation_runner.c, crossval_runner.c and tests/test.c
 * all pass D = days - 1, while simulate_days() integrates t from 0 to
 * days * 24 h — an elapsed time of `days`, not `days - 1`.  The reported RGR
 * is therefore inflated by the factor days / (days - 1): 3.45 % at days = 30.
 * This tool reports both the as-implemented value and the elapsed-correct one
 * so the size of that convention can be seen rather than argued about.
 *
 * Interval RGR
 * ------------
 * For a window [t0, t_end] in days,
 *
 *     RGR(t0, t_end) = ( log B(t_end) - log B(t0) ) / (t_end - t0)
 *
 * with B = B_l + B_r sampled from the same trajectory, so the spin-up variants
 * cost nothing beyond the one run that already produced the horizon points.
 *
 * Integration path
 * ----------------
 * The runs here reproduce simulate_days() exactly — same ODE45 driver, same
 * tolerance, same soil rescaling of lambda_sb, same night-start starch
 * latching, same step budget — with the GUI, the CSV writer and the solver
 * logger switched off and one extra observer added that records the full
 * biomass trajectory plus the individual terms of the growth equations.  It
 * does not call simulate_days() itself because SimulationResult carries only
 * total biomass (result->total_biomass), whereas the growth-term decomposition in
 * part (3) needs leaf and root separately.
 *
 * Outputs (results/horizon/)
 * --------------------------
 *   horizon_sweep.csv        RGR vs total horizon, per photoperiod
 *   spinup_sweep.csv         RGR vs spin-up start, for t_end = 30 and 35 d
 *   conventions.csv          A / B / B* / C: RGR, deviation, recalibrated
 *                            lambda_sb
 *   growth_trace.csv         d(ln B)/dt and the magnitude of every term of
 *                            dB_l/dt over time, per photoperiod
 *   harness_compare.csv      validation_runner vs tests/test.c lambda_sb sets
 */

#define _POSIX_C_SOURCE 200809L

#include "config/config.h"
#include "model/input.h"
#include "model/model.h"
#include "model/state.h"
#include "simulation/dynamics.h"
#include "simulation/observers.h"
#include "simulation/simulation.h"
#include "util/fpclass.h"

#include <math.h>
#include <solver.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#ifndef max_align_t
#define max_align_t long double
#define HORIZON_MAX_ALIGN_T_SHIM
#endif
#include <framework/default_runtime.h>
#include <framework/model_buffers.h>
#include <framework/model_runner.h>
#ifdef HORIZON_MAX_ALIGN_T_SHIM
#undef max_align_t
#undef HORIZON_MAX_ALIGN_T_SHIM
#endif

/* ------------------------------------------------------------------ scenarios */

typedef struct {
    int    photoperiod_h;
    double lit_rgr;      /* Sulpice et al. (2014) Suppl. Table 1, RGR(FW) */
    double lambda_sni;
    double lambda_sb_val;  /* tools/validation_runner.c SCENARIOS[] */
    double lambda_sb_test; /* tests/test.c TEST_RGR macro            */
    double feedback;
} Scenario;

static const Scenario SCENARIOS[] = {
    {  4, 0.068008, 0.16, 0.003843, 0.00344, 0.82 },
    {  6, 0.113516, 0.15, 0.004433, 0.00413, 0.79 },
    {  8, 0.170809, 0.13, 0.005246, 0.00515, 0.67 },
    { 12, 0.260002, 0.08, 0.005850, 0.00578, 0.62 },
};
#define N_SCEN ((int)(sizeof(SCENARIOS) / sizeof(SCENARIOS[0])))

/* Dusk starch targets, Sulpice Suppl. Table S4, read at ZT = photoperiod. */
static const double DUSK_STARCH[N_SCEN] = { 29.6749, 54.2372, 56.8778, 58.7337 };

#define MAX_DAYS 60

/* ------------------------------------------------------------------ trajectory observer */

typedef struct {
    double t;
    double leaf, root, total;
    double starch, sucrose, photosynthesis;
    /* terms of dB_l/dt and dB_r/dt, signed as they enter the derivative */
    double leaf_growth_term, leaf_death_term, leaf_comp_term;
    double root_growth_term, root_death_term, root_comp_term;
    int    gate_active;   /* B_l < min_leaf_biomass, so photosynthesis is forced to 0 */
} TrajPoint;

typedef struct {
    double     next_t;
    double     dt_sample;
    int        step;
    int        max_steps;
    PlantCtx*  rhs_ctx;
    TrajPoint* pts;
} TrajCtx;

static void traj_fill(TrajPoint* p, double t, const Input* s)
{
    const double lsb = (double)s->growth.lambda_sb;
    const double phi_r = (double)s->growth.sucrose_root_allocation;
    const double load = (double)s->carbohydrates.sucrose_loading_frequency;
    const double en = (double)s->carbohydrates.night_efficiency_starch;
    const double bl = (double)s->growth.leaf_biomass;
    const double br = (double)s->growth.root_biomass;

    p->t = t;
    p->leaf = bl;
    p->root = br;
    p->total = bl + br;
    p->starch = (double)s->carbohydrates.starch;
    p->sucrose = (double)s->carbohydrates.sucrose;
    p->photosynthesis = (double)s->photo.photosynthesis;

    p->leaf_growth_term = lsb * (1.0 - phi_r) * load * en * bl;
    p->leaf_death_term = -(double)s->growth.leaf_deathrate * bl;
    p->leaf_comp_term = -(double)s->growth.leaf_competitive_rate * bl * bl;

    p->root_growth_term = lsb * phi_r * load * en * bl;
    p->root_death_term = -(double)s->growth.root_deathrate * br;
    p->root_comp_term = -(double)s->growth.root_competitive_rate * br * br;

    p->gate_active = (bl < (double)s->photo.min_leaf_biomass);
}

static void traj_observer(real_t t, const real_t* x, size_t n, void* vctx)
{
    (void)n;
    TrajCtx* o = (TrajCtx*)vctx;

    if ((double)t + 1e-15 < o->next_t)
        return;
    if (o->step >= o->max_steps)
        return;

    while (o->step < o->max_steps && o->next_t <= (double)t + 1e-15) {
        Input tmp = *(o->rhs_ctx->base);
        simulation_unpack_state(&tmp, x);
        simulation_compute_algebraic(&tmp, (real_t)fmod(o->next_t, 24.0),
                                     o->rhs_ctx->starch_night_start);
        traj_fill(&o->pts[o->step], o->next_t, &tmp);
        o->step++;
        o->next_t += o->dt_sample;
    }
}

/* ------------------------------------------------------------------ step budget */

typedef struct {
    unsigned long steps;
    unsigned long budget;
    int           exhausted;
} BudgetCtx;

static SolverPluginResult budget_on_attempt(void* vctx, const SolverEvent* event)
{
    (void)event;
    BudgetCtx* b = (BudgetCtx*)vctx;
    if (!b)
        return SOLVER_PLUGIN_CONTINUE;
    if (++b->steps >= b->budget) {
        b->exhausted = 1;
        return SOLVER_PLUGIN_STOP;
    }
    return SOLVER_PLUGIN_CONTINUE;
}

/* ------------------------------------------------------------------ one run */

typedef struct {
    TrajPoint* pts;      /* n_pts entries, sample i at t = i * dt hours   */
    int        n_pts;    /* days*24/dt + 1: the final state is appended   */
    double     dt;
    int        days;
    double     fw_start;
    double     fw_end;
    int        diverged;
} Run;

static void run_free(Run* r)
{
    free(r->pts);
    r->pts = NULL;
    r->n_pts = 0;
}

/*
 * Integrate the production model for `days` days and record the trajectory.
 * Mirrors simulate_days() with GUI/CSV/logger disabled; see file header.
 */
static int run_model(const Input* base, const Scenario* sc, double lambda_sb,
                     int days, double dt, double tol, Run* out)
{
    memset(out, 0, sizeof(*out));

    Input input = *base;
    input.photo.photoperiod = (real_t)sc->photoperiod_h;
    input.carbohydrates.lambda_sni = (real_t)sc->lambda_sni;
    input.growth.lambda_sb = (real_t)lambda_sb;
    input.photo.feedback_on_photosynthesis = (real_t)sc->feedback;

    out->fw_start = (double)input.growth.total_biomass;
    out->dt = dt;
    out->days = days;

    /* Soil-dependent rescaling, exactly as simulate_days() applies it. */
    input.growth.lambda_sb = lambda_sb_f(input.growth.lambda_sb,
                                         input.nutrients.nitrogen_soil_content,
                                         input.nutrients.phosphorus_soil_content);

    const int n_samples = (int)llround((double)days * 24.0 / dt);
    out->n_pts = n_samples + 1;
    out->pts = (TrajPoint*)calloc((size_t)out->n_pts, sizeof(TrajPoint));
    if (!out->pts)
        return 1;

    PlantCtx rhs_ctx = {
        .base = &input,
        .starch_night_start = input.carbohydrates.starch,
        .was_light = 0
    };
    {
        Input probe = input;
        update_light_conditions(&probe, REAL(0.0));
        rhs_ctx.was_light = (probe.photo.light != REAL(0.0));
    }

    NanGuardCtx nan_ctx = { 0, 0.0, 0, 0, 0 };
    NightStartObserverCtx night_ctx = { .plant_ctx = &rhs_ctx, .starch_index = X_STARCH };
    TrajCtx traj_ctx = {
        .next_t = 0.0, .dt_sample = dt, .step = 0,
        .max_steps = n_samples, .rhs_ctx = &rhs_ctx, .pts = out->pts
    };

    VectorObserver obs_list[3] = { nan_guard_observer, night_start_observer, traj_observer };
    void* ctx_list[3] = { &nan_ctx, &night_ctx, &traj_ctx };
    ObserverChain chain = { .obs = obs_list, .ctx = ctx_list, .count = 3 };

    SolverStats stats;
    solver_stats_init(&stats);
    SolverPluginManager plugins;
    solver_plugin_manager_init(&plugins);
    SolverDefaultRuntimeConfig rt_cfg = {
        .observer = observer_chain_cb,
        .observer_ctx = &chain,
        .stats = &stats,
        .logger = NULL,
        .async_logger = NULL
    };
    if (solver_default_runtime_init(&plugins, &rt_cfg) != 0)
        plugins.count = 0;

    BudgetCtx budget = { 0, 10000000UL, 0 };
    SolverPlugin budget_plugin = {
        .name = "horizon_step_budget",
        .plugin_ctx = &budget,
        .on_run_start = NULL,
        .on_step_attempt = budget_on_attempt,
        .on_step_accepted = NULL,
        .on_step_rejected = NULL,
        .on_run_finish = NULL
    };
    solver_plugin_manager_add(&plugins, &budget_plugin);

    ModelInterface model = *simulation_plant_model_interface();
    model.model_ctx = &rhs_ctx;
    ModelBuffers buffers = { 0 };
    if (model_buffers_init(&buffers, model.descriptor) != MODEL_STATUS_OK) {
        run_free(out);
        return 1;
    }
    simulation_pack_state(buffers.state, &input);

    SolverRunConfig run_cfg = {
        .type = SOLVER_ODE45,
        .t0 = REAL(0.0),
        .t_end = (real_t)((double)days * 24.0),
        .h_init = (real_t)dt,
        .tol = (real_t)tol,
        .plugins = &plugins
    };
    ModelStatus st = model_run(&model, &buffers, &run_cfg);
    if (st == MODEL_STATUS_OK)
        simulation_unpack_state(&input, buffers.state);
    model_buffers_free(&buffers);

    /* Append the terminal state at t = days*24, which the sampler never
     * reaches (it emits at 0, dt, ..., days*24 - dt). */
    simulation_compute_algebraic(&input, REAL(0.0), rhs_ctx.starch_night_start);
    traj_fill(&out->pts[n_samples], (double)days * 24.0, &input);

    out->fw_end = (double)input.growth.total_biomass;
    out->diverged = (st != MODEL_STATUS_OK) || nan_ctx.tripped || budget.exhausted
                    || pm_isnonfinite(out->fw_end) || (out->fw_end <= 0.0);
    return 0;
}

/* Total biomass at time t_days, read off the sample grid. */
static double biomass_at(const Run* r, double t_days)
{
    long idx = lround(t_days * 24.0 / r->dt);
    if (idx < 0)
        idx = 0;
    if (idx > (long)r->n_pts - 1)
        idx = (long)r->n_pts - 1;
    return r->pts[idx].total;
}

/* RGR over [t0, t_end] in days.  Passing divisor <= 0 uses the elapsed time. */
static double rgr_window(const Run* r, double t0, double t_end, double divisor)
{
    const double b0 = biomass_at(r, t0);
    const double b1 = biomass_at(r, t_end);
    const double d = (divisor > 0.0) ? divisor : (t_end - t0);
    if (b0 <= 0.0 || b1 <= 0.0 || d <= 0.0)
        return NAN;
    return (log(b1) - log(b0)) / d;
}

/* End-of-day starch on the final day: t = (days-1)*24 + photoperiod. */
static double eod_starch(const Run* r, int photoperiod_h)
{
    const double t = (double)(r->days - 1) * 24.0 + (double)photoperiod_h;
    long idx = lround(t / r->dt);
    if (idx < 0)
        idx = 0;
    if (idx > (long)r->n_pts - 1)
        idx = (long)r->n_pts - 1;
    return r->pts[idx].starch;
}

/* ------------------------------------------------------------------ conventions */

typedef struct {
    const char* name;
    int    days;
    double t0;
    double t_end;
    double divisor;    /* <= 0 means use elapsed (t_end - t0) */
    const char* note;
} Convention;

static const Convention CONVENTIONS[] = {
    { "A",  30,  5.0, 30.0, -1.0, "30 d, RGR over days 5..30, elapsed divisor 25" },
    { "B",  30,  0.0, 30.0, 29.0, "30 d, RGR over 0..30, divisor 29 (as implemented today)" },
    { "B*", 30,  0.0, 30.0, -1.0, "30 d, RGR over 0..30, elapsed divisor 30" },
    { "C",  35,  5.0, 35.0, -1.0, "35 d, RGR over days 5..35, elapsed divisor 30" },
};
#define N_CONV ((int)(sizeof(CONVENTIONS) / sizeof(CONVENTIONS[0])))

/* ------------------------------------------------------------------ bisection */

typedef struct {
    double lambda_sb;
    double rgr;
    int    iterations;
    int    status;   /* 0 ok, 1 no bracket, 2 divergent run */
} Calib;

static int rgr_for(const Input* base, const Scenario* sc, const Convention* cv,
                   double lambda_sb, double dt, double tol, double* out_rgr)
{
    Run r;
    if (run_model(base, sc, lambda_sb, cv->days, dt, tol, &r) != 0)
        return 1;
    if (r.diverged) {
        run_free(&r);
        return 2;
    }
    *out_rgr = rgr_window(&r, cv->t0, cv->t_end, cv->divisor);
    run_free(&r);
    return pm_isfinite(*out_rgr) ? 0 : 2;
}

/*
 * Bisect lambda_sb against the published RGR under one convention.
 *
 * Runtime of a single simulation grows with lambda_sb (the biomass, and with
 * it the stiffness of the gas-exchange solve, grows exponentially in it), so
 * the bracket expansion refuses to evaluate lambda_sb > LAMBDA_SB_CAP rather
 * than discovering the cost at run time.  Every calibration in this project
 * lands near 5e-3, three orders of magnitude below the cap.
 */
#define LAMBDA_SB_CAP 1.0

static Calib calibrate(const Input* base, const Scenario* sc, const Convention* cv,
                       double dt, double tol)
{
    Calib c = { NAN, NAN, 0, 1 };
    const double target = sc->lit_rgr;

    double lo = 1.0e-4, hi = 1.0e-2;
    double f_lo, f_hi;

    if (rgr_for(base, sc, cv, lo, dt, tol, &f_lo) != 0)
        return c;
    if (rgr_for(base, sc, cv, hi, dt, tol, &f_hi) != 0)
        return c;

    int expand = 0;
    while (f_hi < target && expand < 12) {
        lo = hi;
        f_lo = f_hi;
        hi *= 2.0;
        if (hi > LAMBDA_SB_CAP)
            return c;
        if (rgr_for(base, sc, cv, hi, dt, tol, &f_hi) != 0)
            return c;
        expand++;
    }
    while (f_lo > target && expand < 24) {
        hi = lo;
        f_hi = f_lo;
        lo *= 0.5;
        if (lo < 1.0e-8)
            return c;
        if (rgr_for(base, sc, cv, lo, dt, tol, &f_lo) != 0)
            return c;
        expand++;
    }
    if (!(f_lo <= target && target <= f_hi))
        return c;

    double mid = 0.5 * (lo + hi), f_mid = NAN;
    for (int it = 0; it < 60; ++it) {
        mid = 0.5 * (lo + hi);
        if (rgr_for(base, sc, cv, mid, dt, tol, &f_mid) != 0) {
            c.status = 2;
            return c;
        }
        c.iterations = it + 1;
        if (f_mid < target)
            lo = mid;
        else
            hi = mid;
        if (fabs(f_mid - target) < 1.0e-10 || (hi - lo) < 1.0e-12)
            break;
    }
    c.lambda_sb = mid;
    c.rgr = f_mid;
    c.status = 0;
    return c;
}

/* ------------------------------------------------------------------ main */

static void ensure_dir(const char* p) { mkdir(p, 0775); }

static void usage(const char* a0)
{
    printf("Usage: %s [--config path] [--dt H] [--tol X] [--out-dir DIR]\n", a0);
}

int main(int argc, char** argv)
{
    const char* config_path = "config/default.config";
    const char* out_dir = "results/horizon";
    double dt = 0.1;
    double tol = 1e-8;

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--config") == 0 && i + 1 < argc)
            config_path = argv[++i];
        else if (strcmp(argv[i], "--dt") == 0 && i + 1 < argc)
            dt = atof(argv[++i]);
        else if (strcmp(argv[i], "--tol") == 0 && i + 1 < argc)
            tol = atof(argv[++i]);
        else if (strcmp(argv[i], "--out-dir") == 0 && i + 1 < argc)
            out_dir = argv[++i];
        else {
            usage(argv[0]);
            return (strcmp(argv[i], "--help") == 0) ? 0 : 1;
        }
    }

    /* config_load_file only assigns keys it finds, so both structures must be
     * fully initialised first or unlisted fields silently stay zero. */
    SimulationConfig cfg = simulation_default_config();
    Input base_input = generate_input();
    if (config_load_file(config_path, &cfg, &base_input) != 0) {
        fprintf(stderr, "Failed to load config: %s\n", config_path);
        return 1;
    }

    ensure_dir("results");
    ensure_dir(out_dir);
    char path[512];

    /* ---------------------------------------------------------- (a) horizon sweep */

    static const int HORIZONS[] = { 15, 20, 25, 30, 35, 40, 50 };
    const int n_hor = (int)(sizeof(HORIZONS) / sizeof(HORIZONS[0]));

    snprintf(path, sizeof(path), "%s/horizon_sweep.csv", out_dir);
    FILE* fh = fopen(path, "w");
    if (!fh) { fprintf(stderr, "cannot write %s\n", path); return 1; }
    fprintf(fh, "photoperiod_h,days,lambda_sb,rgr_impl_div_days_minus_1,"
                "rgr_elapsed_div_days,lit_rgr,rel_dev_impl_pct,rel_dev_elapsed_pct,"
                "biomass_start,biomass_end,leaf_end,root_end,eod_starch_last_day,"
                "dusk_starch_target,diverged\n");

    for (int s = 0; s < N_SCEN; ++s) {
        const Scenario* sc = &SCENARIOS[s];
        for (int h = 0; h < n_hor; ++h) {
            const int days = HORIZONS[h];
            Run r;
            if (run_model(&base_input, sc, sc->lambda_sb_val, days, dt, tol, &r) != 0)
                continue;
            const double rgr_impl = rgr_window(&r, 0.0, (double)days, (double)(days - 1));
            const double rgr_el = rgr_window(&r, 0.0, (double)days, -1.0);
            const int last = r.n_pts - 1;
            fprintf(fh, "%d,%d,%.8g,%.8f,%.8f,%.6f,%.4f,%.4f,"
                        "%.8g,%.8g,%.8g,%.8g,%.6f,%.4f,%d\n",
                    sc->photoperiod_h, days, sc->lambda_sb_val,
                    rgr_impl, rgr_el, sc->lit_rgr,
                    100.0 * (rgr_impl - sc->lit_rgr) / sc->lit_rgr,
                    100.0 * (rgr_el - sc->lit_rgr) / sc->lit_rgr,
                    r.pts[0].total, r.pts[last].total,
                    r.pts[last].leaf, r.pts[last].root,
                    eod_starch(&r, sc->photoperiod_h), DUSK_STARCH[s], r.diverged);
            run_free(&r);
        }
        printf("horizon sweep done: %d h\n", sc->photoperiod_h);
        fflush(stdout);
    }
    fclose(fh);

    /* ---------------------------------------------------------- (b) spin-up sweep */

    static const double T0S[] = { 0.0, 2.0, 5.0, 7.0, 10.0 };
    const int n_t0 = (int)(sizeof(T0S) / sizeof(T0S[0]));
    static const int T_ENDS[] = { 30, 35 };

    snprintf(path, sizeof(path), "%s/spinup_sweep.csv", out_dir);
    FILE* fs = fopen(path, "w");
    if (!fs) { fprintf(stderr, "cannot write %s\n", path); return 1; }
    fprintf(fs, "photoperiod_h,t_end_days,t0_days,window_days,lambda_sb,rgr,lit_rgr,"
                "rel_dev_pct,biomass_t0,biomass_t_end,leaf_t_end,root_t_end,"
                "eod_starch_last_day,dusk_starch_target,diverged\n");

    for (int s = 0; s < N_SCEN; ++s) {
        const Scenario* sc = &SCENARIOS[s];
        for (int e = 0; e < 2; ++e) {
            const int days = T_ENDS[e];
            Run r;
            if (run_model(&base_input, sc, sc->lambda_sb_val, days, dt, tol, &r) != 0)
                continue;
            for (int k = 0; k < n_t0; ++k) {
                const double t0 = T0S[k];
                const double rgr = rgr_window(&r, t0, (double)days, -1.0);
                const int last = r.n_pts - 1;
                fprintf(fs, "%d,%d,%.1f,%.1f,%.8g,%.8f,%.6f,%.4f,"
                            "%.8g,%.8g,%.8g,%.8g,%.6f,%.4f,%d\n",
                        sc->photoperiod_h, days, t0, (double)days - t0,
                        sc->lambda_sb_val, rgr, sc->lit_rgr,
                        100.0 * (rgr - sc->lit_rgr) / sc->lit_rgr,
                        biomass_at(&r, t0), r.pts[last].total,
                        r.pts[last].leaf, r.pts[last].root,
                        eod_starch(&r, sc->photoperiod_h), DUSK_STARCH[s], r.diverged);
            }
            run_free(&r);
        }
        printf("spin-up sweep done: %d h\n", sc->photoperiod_h);
        fflush(stdout);
    }
    fclose(fs);

    /* ---------------------------------------------------------- (c) growth trace */

    snprintf(path, sizeof(path), "%s/growth_trace.csv", out_dir);
    FILE* fg = fopen(path, "w");
    if (!fg) { fprintf(stderr, "cannot write %s\n", path); return 1; }
    fprintf(fg, "photoperiod_h,t_days,leaf,root,total,dlnB_dt_per_day,"
                "leaf_growth_term,leaf_death_term,leaf_comp_term,"
                "comp_over_growth,death_over_growth,"
                "root_growth_term,root_death_term,root_comp_term,"
                "starch,sucrose,photosynthesis,gate_active\n");

    for (int s = 0; s < N_SCEN; ++s) {
        const Scenario* sc = &SCENARIOS[s];
        Run r;
        if (run_model(&base_input, sc, sc->lambda_sb_val, 35, dt, tol, &r) != 0)
            continue;
        /* Instantaneous d(ln B)/dt from the recorded derivative terms, so it is
         * the model's own rate and not a finite difference of the samples. */
        const int stride = (int)llround(1.0 / dt);   /* one point per hour */
        for (int i = 0; i < r.n_pts; i += (stride > 0 ? stride : 1)) {
            const TrajPoint* p = &r.pts[i];
            const double dB = p->leaf_growth_term + p->leaf_death_term + p->leaf_comp_term
                            + p->root_growth_term + p->root_death_term + p->root_comp_term;
            const double dlnB = (p->total > 0.0) ? 24.0 * dB / p->total : NAN;
            const double g = p->leaf_growth_term;
            fprintf(fg, "%d,%.6f,%.10g,%.10g,%.10g,%.10g,"
                        "%.10g,%.10g,%.10g,%.6e,%.6e,"
                        "%.10g,%.10g,%.10g,%.8g,%.8g,%.8g,%d\n",
                    sc->photoperiod_h, p->t / 24.0, p->leaf, p->root, p->total, dlnB,
                    p->leaf_growth_term, p->leaf_death_term, p->leaf_comp_term,
                    (g > 0.0) ? fabs(p->leaf_comp_term) / g : 0.0,
                    (g > 0.0) ? fabs(p->leaf_death_term) / g : 0.0,
                    p->root_growth_term, p->root_death_term, p->root_comp_term,
                    p->starch, p->sucrose, p->photosynthesis, p->gate_active);
        }
        run_free(&r);
        printf("growth trace done: %d h\n", sc->photoperiod_h);
        fflush(stdout);
    }
    fclose(fg);

    /* ---------------------------------------------------------- (d) conventions */

    snprintf(path, sizeof(path), "%s/conventions.csv", out_dir);
    FILE* fc = fopen(path, "w");
    if (!fc) { fprintf(stderr, "cannot write %s\n", path); return 1; }
    fprintf(fc, "convention,note,photoperiod_h,days,t0_days,t_end_days,divisor,"
                "lambda_sb_frozen,rgr_frozen,lit_rgr,rel_dev_pct,"
                "lambda_sb_recalibrated,rgr_recalibrated,recal_over_frozen,"
                "calib_status,calib_iterations\n");

    for (int c = 0; c < N_CONV; ++c) {
        const Convention* cv = &CONVENTIONS[c];
        for (int s = 0; s < N_SCEN; ++s) {
            const Scenario* sc = &SCENARIOS[s];
            double rgr_frozen = NAN;
            rgr_for(&base_input, sc, cv, sc->lambda_sb_val, dt, tol, &rgr_frozen);
            Calib cal = calibrate(&base_input, sc, cv, dt, tol);
            fprintf(fc, "%s,\"%s\",%d,%d,%.1f,%.1f,%.4f,"
                        "%.8g,%.8f,%.6f,%.4f,%.8g,%.8f,%.6f,%d,%d\n",
                    cv->name, cv->note, sc->photoperiod_h, cv->days, cv->t0, cv->t_end,
                    (cv->divisor > 0.0) ? cv->divisor : (cv->t_end - cv->t0),
                    sc->lambda_sb_val, rgr_frozen, sc->lit_rgr,
                    100.0 * (rgr_frozen - sc->lit_rgr) / sc->lit_rgr,
                    cal.lambda_sb, cal.rgr,
                    cal.lambda_sb / sc->lambda_sb_val,
                    cal.status, cal.iterations);
            fflush(fc);
        }
        printf("convention %s done\n", cv->name);
        fflush(stdout);
    }
    fclose(fc);

    /* ---------------------------------------------------------- (e) harness compare */

    snprintf(path, sizeof(path), "%s/harness_compare.csv", out_dir);
    FILE* fx = fopen(path, "w");
    if (!fx) { fprintf(stderr, "cannot write %s\n", path); return 1; }
    fprintf(fx, "photoperiod_h,source,lambda_sb,days,rgr_impl_div29,rgr_elapsed_div30,"
                "lit_rgr,rel_dev_impl_pct,rel_dev_elapsed_pct,"
                "passes_10pct_impl,passes_10pct_elapsed,"
                "biomass_end,eod_starch_last_day,dusk_starch_target\n");

    const char* SRC_NAME[2] = { "validation_runner", "tests_test_c" };
    for (int s = 0; s < N_SCEN; ++s) {
        const Scenario* sc = &SCENARIOS[s];
        const double lsbs[2] = { sc->lambda_sb_val, sc->lambda_sb_test };
        for (int k = 0; k < 2; ++k) {
            Run r;
            if (run_model(&base_input, sc, lsbs[k], 30, dt, tol, &r) != 0)
                continue;
            const double rgr_impl = rgr_window(&r, 0.0, 30.0, 29.0);
            const double rgr_el = rgr_window(&r, 0.0, 30.0, -1.0);
            fprintf(fx, "%d,%s,%.8g,30,%.8f,%.8f,%.6f,%.4f,%.4f,%d,%d,%.8g,%.6f,%.4f\n",
                    sc->photoperiod_h, SRC_NAME[k], lsbs[k],
                    rgr_impl, rgr_el, sc->lit_rgr,
                    100.0 * (rgr_impl - sc->lit_rgr) / sc->lit_rgr,
                    100.0 * (rgr_el - sc->lit_rgr) / sc->lit_rgr,
                    (fabs(rgr_impl - sc->lit_rgr) < sc->lit_rgr / 10.0) ? 1 : 0,
                    (fabs(rgr_el - sc->lit_rgr) < sc->lit_rgr / 10.0) ? 1 : 0,
                    r.fw_end, eod_starch(&r, sc->photoperiod_h), DUSK_STARCH[s]);
            run_free(&r);
        }
        printf("harness compare done: %d h\n", sc->photoperiod_h);
        fflush(stdout);
    }
    fclose(fx);

    printf("horizon study complete; CSVs in %s\n", out_dir);
    return 0;
}

/*
 * convergence_study.c — numerical convergence of the reported results
 *
 * Purpose
 * -------
 * Reviewer issue #15 asks whether the reported quantities are numerically
 * converged, or whether they are artefacts of the integration step size and of
 * the solver tolerances.  The question is answered by refinement rather than by
 * assertion: every reported quantity is recomputed over a sequence of step
 * sizes and tolerances that brackets the production setting, and the change
 * between successive levels is what is reported.
 *
 * WHAT dt_hours ACTUALLY IS (this governs the whole design of the study)
 * ---------------------------------------------------------------------
 * The production integrator is the adaptive ODE45 backend
 * (config/default.config: simulation.solver_type = ODE45,
 * simulation.solver_tolerance = 1e-8), not fixed-step RK4.  Under an adaptive
 * integrator cfg.dt_hours does two unrelated jobs
 * (src/simulation/simulation.c): it is passed as SolverRunConfig.h_init, i.e.
 * the *initial* trial step that the controller immediately overrides, and it is
 * passed as runtime->sample_ctx.dt_sample, i.e. the interval at which
 * sim_sampling_observer records the trajectory.  It is NOT the integration
 * step.  Sweeping it therefore refines the OUTPUT GRID, not the
 * discretisation, and the accuracy of the production trajectory is controlled
 * by solver_tolerance alone.
 *
 * An earlier version of this study swept dt through a fixed-step RK4 driver
 * local to this file and reported the result as the convergence of the
 * published quantities.  That conflated three different things — the
 * production integrator with a non-production one, the integration step with
 * the sampling interval, and a discretisation error with an output-grid error.
 * The sweeps below separate them explicitly.
 *
 * Sweeps, for each of the four photoperiods (4, 6, 8, 12 h)
 * --------------------------------------------------------
 *   (a) PRODUCTION tolerance refinement.  ODE45 through simulate_days(), dt
 *       held at the production 0.1 h, solver_tolerance swept over
 *       {1e-4 ... 1e-12}.  This is the production-relevant convergence test:
 *       it refines the quantity that actually controls the error.  The
 *       observed order is estimated by Richardson extrapolation on successive
 *       triples,
 *
 *           p = log2( |f(l) - f(l+1)| / |f(l+1) - f(l+2)| ),
 *
 *       which requires no exact solution.  Because successive tolerance levels
 *       are a factor of 100 apart, the same differences are also reported as an
 *       order with respect to the actual refinement ratio,
 *       log_r(d_coarse/d_fine), in convergence_richardson.csv.
 *
 *   (b) PRODUCTION dt refinement at fixed tolerance.  ODE45, tolerance held at
 *       the production 1e-8, dt swept over {0.4 ... 0.00625} h.  Since dt is
 *       h_init plus the sampling interval, whatever this sweep moves is an
 *       output-grid effect (the end-of-day sample is the first accepted step at
 *       or past the sample time, with no interpolation), not a discretisation
 *       effect.  Comparing (a) against (b) at the same production point is what
 *       attributes the error to its true source; the attribution is tabulated
 *       in convergence_error_separation.csv.
 *
 *   (c) RK4 fixed-step refinement — BACKEND VERIFICATION ONLY, NOT the
 *       production configuration.  The fixed-step RK4 driver local to this file
 *       with dt halved over seven levels.  Here dt genuinely is the integration
 *       step, so this sweep measures a discretisation error and a meaningful
 *       observed order, and it is the only path that can carry the exact
 *       carbon-balance bookkeeping (ten flux accumulators advanced by the same
 *       RK4 weights).  It is retained because it is a legitimate verification
 *       of the integrator and the only available conservation check — but no
 *       number from it describes the published runs, which are ODE45.
 *
 *   (d) Ci solver tolerance.  eps in {1e-3, 1e-4, 1e-6, 1e-8, 1e-10} applied to
 *       the legacy fixed-point scheme through iterate_ci_tol(eps_ci, eps_an);
 *       1e-4 is the production value of the historical single-epsilon code.
 *       The default analytic solver is closed-form and has no convergence
 *       tolerance at all, so this sweep necessarily runs in
 *       CI_SOLVER_FIXED_POINT mode.  That the production quantity is
 *       tolerance-free under the default is itself part of the answer.
 *
 *   (e) Ci solver mode.  Analytic (default) against legacy fixed-point at
 *       psi_leaf = -0.6 MPa, which is the value all production runs use.  At
 *       this water potential the fixed-point map is contractive and the two
 *       solvers should agree closely; that agreement is the result worth
 *       reporting, because it establishes that replacing the solver did not
 *       move any published number.  The same comparison is repeated at
 *       psi_leaf = -1.2 MPa, past the non-convergence threshold documented in
 *       docs/revision/math_audit.md section 4.1, to show what the analytic
 *       solver buys outside the production regime.
 *
 *   (f) Integrator identity.  The production config selects the adaptive
 *       ODE45 backend, while several production tools (sensitivity_runner,
 *       oat_sweep, multi_2d_sweep) override to fixed-step RK4.  Both are run
 *       through simulate_days() at production settings so the report can state
 *       whether the choice of backend is visible in the reported quantities.
 *
 *   (g) Backend sanity against an exact solution, dy/dt = -y.
 *
 * Reported quantities
 * -------------------
 * For every configuration: RGR, final total biomass, end-of-day starch,
 * end-of-day cytosolic sucrose, end-of-day total sucrose (cytosol + vacuole),
 * and the relative residual of the carbohydrate carbon balance.
 *
 * The carbon-balance residual is obtained the same way as in
 * tools/carbon_balance.c: ten flux accumulators are appended to the eleven
 * states and advanced by the same RK4 weights and the same step size, so the
 * budget is exact to the arithmetic of the integrator.  Carrying it through the
 * dt sweep is deliberate — it distinguishes a discretisation error in the
 * trajectory (which shrinks with dt) from a conservation defect in the flux
 * bookkeeping (which would not).
 *
 * Why no measured order is expected to be 4
 * ----------------------------------------
 * Under ODE45 there is no "order" in the classical sense to measure at all:
 * sweep (a) refines a tolerance, not a step, so its Richardson exponent is a
 * diagnostic of how the error responds to the requested accuracy and has no
 * formal target value; and sweep (b) refines an output grid, so its exponent
 * describes sampling, not the method.  Only sweep (c) refines an integration
 * step — and even there 4 is not the expectation, for the reason that follows.
 * update_light_conditions() switches the light flag at sunrise and sunset, and
 * starch_degradation() is multiplied by (1 - light).  The right-hand side is
 * therefore discontinuous in t twice per simulated day.  Classical RK4 attains
 * its design order only for a right-hand side that is C^4 across the step; a
 * step whose stages straddle a jump commits an O(h) local error irrespective of
 * the method.  Because the step grid used here divides both 24 h and every
 * photoperiod exactly, no step straddles a transition, but the stage evaluated
 * at exactly t = sunset lands on the dark branch while the step it belongs to
 * is a light-period step.  The result is one O(h^2) local defect per transition
 * rather than O(h^5), and the global order is limited by those transitions and
 * not by the interior of the day.  The measured p is reported as measured; the
 * dt sweep additionally reports whether the residual change between the
 * production step and the next two refinements is small in absolute terms,
 * which is the practically relevant statement regardless of the order.
 *
 * Outputs (results/convergence/)
 * ------------------------------
 *   convergence_solver_tol.csv        (a) ODE45, dt = 0.1 h fixed, tol swept
 *                                     -- THE PRODUCTION-RELEVANT TEST
 *   convergence_dt_ode45.csv          (b) ODE45, tol = 1e-8 fixed, dt swept
 *                                     -- output-grid effect
 *   convergence_error_separation.csv  (a) vs (b) side by side, per photoperiod
 *   convergence_richardson.csv        Richardson orders for (a), (b) and (c),
 *                                     one row per (sweep, photoperiod, metric,
 *                                     triple)
 *   convergence_dt.csv                (c) fixed-step RK4 dt refinement
 *                                     -- VERIFICATION ONLY, NOT PRODUCTION
 *   convergence_dt_order.csv          (c) Richardson orders (retained for
 *                                     backward compatibility; also in
 *                                     convergence_richardson.csv)
 *   convergence_tolerance.csv         (d) one row per (photoperiod, eps_ci)
 *   convergence_solver_mode.csv       (e) one row per (photoperiod, psi, mode)
 *   convergence_integrator.csv        (f) one row per (photoperiod, backend)
 *   convergence_integrator_exact.csv  (g) one row per (backend, tol)
 *
 * Usage
 * -----
 *   ./convergence_study [--config <path>] [--out <dir>] [--days N]
 *
 * Must be run with the repository root as the working directory: the default
 * config path is relative.
 */

#include "config/config.h"
#include "model/input.h"
#include "model/model.h"
#include "model/state.h"
#include "simulation/dynamics.h"
#include "simulation/simulation.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

/* ------------------------------------------------------------------ scenarios */

/* Per-photoperiod parameters, identical to the SCENARIOS table in
 * tools/validation_runner.c so that this study refines the published
 * configuration rather than some neighbouring one. */
typedef struct {
    int    photoperiod_h;
    real_t lambda_sni;
    real_t lambda_sb;
    real_t feedback;      /* lambda_c: feedback_on_photosynthesis */
} Scenario;

static const Scenario SCENARIOS[] = {
    {  4, REAL(0.16), REAL(0.003843), REAL(0.82) },
    {  6, REAL(0.15), REAL(0.004433), REAL(0.79) },
    {  8, REAL(0.13), REAL(0.005246), REAL(0.67) },
    { 12, REAL(0.08), REAL(0.005850), REAL(0.62) },
};
#define NPP ((int)(sizeof(SCENARIOS) / sizeof(SCENARIOS[0])))

/* ------------------------------------------------------------------ flux accumulators */

/* Layout mirrors tools/carbon_balance.c. */
typedef enum {
    F_ASSIM = 0, F_TO_STARCH, F_TO_SUCROSE, F_DEGRADE, F_UPTAKE,
    F_TRANSPORT, F_MAINT, F_LOAD, F_VAC_IN, F_VAC_OUT, F_DIM
} FluxIndex;

#define AUG_DIM (X_DIM + F_DIM)

static void augmented_rhs(real_t t, const real_t* x, real_t* dx, PlantCtx* ctx)
{
    simulation_plant_rhs(t, x, dx, ctx);

    Input s = *(ctx->base);
    simulation_unpack_state(&s, x);
    s.growth.total_biomass = s.growth.leaf_biomass + s.growth.root_biomass;
    simulation_compute_algebraic(&s, fmod((double)t, 24.0), ctx->starch_night_start);

    const real_t A_n = s.photo.photosynthesis;
    const real_t k_s = s.carbohydrates.starch_partition_coeff;
    const real_t suc = s.carbohydrates.sucrose;
    const real_t vac = s.carbohydrates.sucrose_vacuole;

    real_t* f = dx + X_DIM;
    f[F_ASSIM]      = A_n;
    f[F_TO_STARCH]  = k_s * A_n;
    f[F_TO_SUCROSE] = (REAL(1.0) - k_s) * A_n;
    f[F_DEGRADE]    = s.carbohydrates.starch_degradation_rate;
    f[F_UPTAKE]     = s.carbohydrates.uptake_cost;
    f[F_TRANSPORT]  = s.carbohydrates.transport_cost;
    f[F_MAINT]      = s.carbohydrates.respiration_frequency * suc;
    f[F_LOAD]       = s.carbohydrates.sucrose_loading_frequency
                    * s.carbohydrates.night_efficiency_starch;
    f[F_VAC_IN]     = s.carbohydrates.k_vacuole_in  * suc;
    f[F_VAC_OUT]    = s.carbohydrates.k_vacuole_out * vac;
}

/* Exact-solution test problem for sweep (e): dy/dt = -y, y(t) = exp(-t). */
static void linear_decay_rhs(real_t t, const real_t* y, real_t* dy,
                             size_t n, void* ctx)
{
    (void)t; (void)ctx;
    for (size_t i = 0; i < n; ++i) dy[i] = -y[i];
}

static void rk4_step(real_t t, real_t h, real_t* x, PlantCtx* ctx,
                     real_t* k1, real_t* k2, real_t* k3, real_t* k4, real_t* tmp)
{
    augmented_rhs(t, x, k1, ctx);
    for (size_t i = 0; i < AUG_DIM; ++i) tmp[i] = x[i] + REAL(0.5) * h * k1[i];
    augmented_rhs(t + REAL(0.5) * h, tmp, k2, ctx);
    for (size_t i = 0; i < AUG_DIM; ++i) tmp[i] = x[i] + REAL(0.5) * h * k2[i];
    augmented_rhs(t + REAL(0.5) * h, tmp, k3, ctx);
    for (size_t i = 0; i < AUG_DIM; ++i) tmp[i] = x[i] + h * k3[i];
    augmented_rhs(t + h, tmp, k4, ctx);
    for (size_t i = 0; i < AUG_DIM; ++i)
        x[i] += (h / REAL(6.0)) * (k1[i] + REAL(2.0) * k2[i] + REAL(2.0) * k3[i] + k4[i]);
}

/* ------------------------------------------------------------------ metrics */

typedef struct {
    double rgr;
    double final_biomass;
    double eod_starch;
    double eod_sucrose;
    double eod_sucrose_total;
    double carbon_residual_rel;
    double carbon_residual_abs;
    double int_assimilation;
    long   ci_calls;
    long   ci_nonconvergent;
    double ci_mean_iterations;
    long   steps_per_day;
    /* Non-zero if the run was classified as divergent (non-finite state, or the
     * step budget exhausted by a collapsed adaptive step).  Only the
     * simulate_days() path can set this; the local RK4 driver always reports 0.
     * Tight tolerances are exactly where a stall is plausible, so every row of
     * the ODE45 sweeps carries this flag rather than silently reporting a
     * truncated trajectory as a converged one. */
    int    diverged;
} Metrics;

/* The five quantities whose convergence is tracked, in a fixed order so that
 * the order-estimation loop can address them generically. */
#define NMETRIC 5
static const char* METRIC_NAME[NMETRIC] = {
    "rgr", "final_biomass", "eod_starch", "eod_sucrose", "eod_sucrose_total"
};

static double metric_get(const Metrics* m, int k)
{
    switch (k) {
    case 0:  return m->rgr;
    case 1:  return m->final_biomass;
    case 2:  return m->eod_starch;
    case 3:  return m->eod_sucrose;
    default: return m->eod_sucrose_total;
    }
}

static Input scenario_input(const Input* base, const Scenario* sc)
{
    Input in = *base;
    in.photo.photo_mode                 = (int)PHOTO_MODE_FARQUHAR_JARVIS;
    in.photo.photoperiod                = (real_t)sc->photoperiod_h;
    in.carbohydrates.lambda_sni         = sc->lambda_sni;
    in.growth.lambda_sb                 = sc->lambda_sb;
    in.photo.feedback_on_photosynthesis = sc->feedback;
    return in;
}

/*
 * Fixed-step RK4 integration of the augmented system, reproducing the
 * production model exactly: the same soil rescaling of lambda_sb that
 * simulate_days() applies, and the same latching of the night-start starch
 * reference at each light -> dark transition.
 *
 * The end-of-day sample is taken at t = (days-1)*24 + photoperiod.  All step
 * sizes used here divide both 24 and every photoperiod exactly, so that sample
 * lands on a grid point and no interpolation is involved; the assertion is
 * checked at run time and a mismatch is reported rather than silently
 * interpolated.
 */
static int run_rk4(const Input* base, const Scenario* sc, int days, real_t dt,
                   Metrics* out)
{
    Input input = scenario_input(base, sc);
    input.growth.lambda_sb = lambda_sb_f(
        input.growth.lambda_sb,
        input.nutrients.nitrogen_soil_content,
        input.nutrients.phosphorus_soil_content);

    const real_t fw_start = input.growth.leaf_biomass + input.growth.root_biomass;

    PlantCtx ctx = { .base = &input,
                     .starch_night_start = input.carbohydrates.starch,
                     .was_light = 0 };
    {
        Input probe = input;
        update_light_conditions(&probe, REAL(0.0));
        ctx.was_light = (probe.photo.light != REAL(0.0));
    }

    real_t x[AUG_DIM];
    memset(x, 0, sizeof(x));
    {
        PlantState st = plant_state_from_input(&input);
        plant_state_pack(x, &st);
    }
    const double a0 = (double)x[X_STARCH];
    const double s0 = (double)x[X_SUCROSE];
    const double v0 = (double)x[X_SUCROSE_VAC];

    const long steps_day = (long)llround(24.0 / (double)dt);
    if (fabs((double)steps_day * (double)dt - 24.0) > 1e-12) {
        fprintf(stderr, "dt = %g does not divide 24 h exactly\n", (double)dt);
        return 1;
    }
    const long eod_step = (long)(days - 1) * steps_day
                        + (long)llround((double)sc->photoperiod_h / (double)dt);
    if (fabs((double)llround((double)sc->photoperiod_h / (double)dt) * (double)dt
             - (double)sc->photoperiod_h) > 1e-12) {
        fprintf(stderr, "dt = %g does not divide photoperiod %d exactly\n",
                (double)dt, sc->photoperiod_h);
        return 1;
    }

    int    eod_seen   = 0;
    double eod_starch = 0.0, eod_sucrose = 0.0, eod_vac = 0.0;

    real_t k1[AUG_DIM], k2[AUG_DIM], k3[AUG_DIM], k4[AUG_DIM], tmp[AUG_DIM];

    ci_solver_reset_diagnostics();

    long step = 0;
    for (int d = 0; d < days; ++d) {
        for (long k = 0; k < steps_day; ++k, ++step) {
            const real_t t = (real_t)((double)d * 24.0 + (double)k * (double)dt);

            /* End-of-day sample: state at exactly t = (days-1)*24 + P, i.e.
             * the instant the light goes off on the final day. */
            if (step == eod_step) {
                eod_starch  = (double)x[X_STARCH];
                eod_sucrose = (double)x[X_SUCROSE];
                eod_vac     = (double)x[X_SUCROSE_VAC];
                eod_seen    = 1;
            }

            Input probe = input;
            update_light_conditions(&probe, (real_t)fmod((double)t, 24.0));
            const int is_light = (probe.photo.light != REAL(0.0));
            if (ctx.was_light && !is_light)
                ctx.starch_night_start = x[X_STARCH];
            ctx.was_light = is_light;

            rk4_step(t, dt, x, &ctx, k1, k2, k3, k4, tmp);
        }
    }

    if (!eod_seen) {
        fprintf(stderr, "end-of-day sample never reached (dt = %g)\n", (double)dt);
        return 1;
    }

    const CiSolverDiagnostics diag = ci_solver_get_diagnostics();

    const double a1 = (double)x[X_STARCH];
    const double s1 = (double)x[X_SUCROSE];
    const double v1 = (double)x[X_SUCROSE_VAC];
    const double fw_end = (double)(x[X_LEAF_BIOMASS] + x[X_ROOT_BIOMASS]);

    /* Carbohydrate budget over [0, T]: pool change minus net flux. */
    const double pool_change = (a1 + s1 + v1) - (a0 + s0 + v0);
    const double net = (double)x[X_DIM + F_ASSIM]
                     - ((double)x[X_DIM + F_UPTAKE]
                      + (double)x[X_DIM + F_TRANSPORT]
                      + (double)x[X_DIM + F_MAINT]
                      + (double)x[X_DIM + F_LOAD]);
    const double res = pool_change - net;

    out->rgr = log(fw_end / (double)fw_start) / (double)(days - 1);
    out->final_biomass       = fw_end;
    out->eod_starch          = eod_starch;
    out->eod_sucrose         = eod_sucrose;
    out->eod_sucrose_total   = eod_sucrose + eod_vac;
    out->int_assimilation    = (double)x[X_DIM + F_ASSIM];
    out->carbon_residual_abs = res;
    out->carbon_residual_rel = (out->int_assimilation > 0.0)
                             ? res / out->int_assimilation : 0.0;
    out->ci_calls            = diag.calls;
    out->ci_nonconvergent    = diag.nonconvergent;
    out->ci_mean_iterations  = (diag.calls > 0)
                             ? (double)diag.iterations_total / (double)diag.calls : 0.0;
    out->steps_per_day       = steps_day;
    out->diverged            = 0;   /* fixed-step driver cannot stall */
    return 0;
}

/*
 * Production path: the same quantities obtained through simulate_days(), so the
 * comparison between integrator backends (sweep d) exercises the code the
 * published runs actually used rather than this tool's own RK4 driver.  The
 * carbon-balance residual is not available on this path — simulate_days()
 * integrates the eleven states only — and is reported as not-applicable.
 */
static int run_production(const SimulationConfig* cfg_in, const Input* base,
                          const Scenario* sc, Metrics* out)
{
    Input in = scenario_input(base, sc);
    SimulationConfig cfg = *cfg_in;
    cfg.gui_enabled = 0;
    cfg.write_csv   = 0;

    const real_t fw_start = in.growth.total_biomass;

    ci_solver_reset_diagnostics();

    SimulationResult r;
    memset(&r, 0, sizeof(r));
    simulate_days(&cfg, &in, &r);

    const CiSolverDiagnostics diag = ci_solver_get_diagnostics();

    const long steps_day  = (long)llround(24.0 / (double)cfg.dt_hours);
    long last_start = (long)r.filled_steps - steps_day;
    if (last_start < 0) last_start = 0;
    long eod = last_start + (long)llround((double)sc->photoperiod_h / (double)cfg.dt_hours);
    if (eod >= (long)r.filled_steps) eod = (long)r.filled_steps - 1;

    out->rgr = (double)compute_rgr(fw_start, in.growth.total_biomass,
                                   (real_t)(cfg.days - 1));
    out->final_biomass     = (double)in.growth.total_biomass;
    out->eod_starch        = (r.filled_steps > 0) ? (double)r.starch[eod] : 0.0;
    out->eod_sucrose       = (r.filled_steps > 0) ? (double)r.sucrose[eod] : 0.0;
    out->eod_sucrose_total = (r.filled_steps > 0) ? (double)r.sucrose_total[eod] : 0.0;
    out->int_assimilation    = 0.0;
    out->carbon_residual_abs = NAN;
    out->carbon_residual_rel = NAN;
    out->ci_calls            = diag.calls;
    out->ci_nonconvergent    = diag.nonconvergent;
    out->ci_mean_iterations  = (diag.calls > 0)
                             ? (double)diag.iterations_total / (double)diag.calls : 0.0;
    out->steps_per_day       = steps_day;
    out->diverged            = r.diverged;

    simulation_result_free(&r);
    return 0;
}

/* ------------------------------------------------------------------ helpers */

static void ensure_dir(const char* path) { mkdir(path, 0775); }

static double rel_dev(double a, double b)
{
    const double d = fabs(b);
    return (d > 0.0) ? fabs(a - b) / d : fabs(a - b);
}

#define NDT 7

/* ODE45 sweep (a): tolerance levels at the production dt.  1e-8 is production;
 * the sweep brackets it four orders below and four above. */
static const double STOL[] = { 1e-4, 1e-6, 1e-8, 1e-10, 1e-12 };
#define NSTOL ((int)(sizeof(STOL) / sizeof(STOL[0])))
#define STOL_PROD_IDX 2
/* Successive tolerance levels differ by a factor of 100, so an order estimated
 * from a triple must use log_100, not log_2. */
#define STOL_RATIO 100.0

/* ODE45 sweep (b): dt levels at the production tolerance.  Same grid as the RK4
 * sweep (c) so that the two are directly comparable level for level. */
#define NDTP 7

/*
 * Richardson order from a triple of successive refinement levels:
 *
 *     p = log_r( |f_0 - f_1| / |f_1 - f_2| ),
 *
 * where r is the refinement ratio between levels (2 for a halved step, 100 for
 * a tolerance tightened by two decades).  NAN when either difference vanishes:
 * with three identical values there is no order to extract, and reporting 0 or
 * a large number would both be inventions.
 */
static double richardson_order(double f0, double f1, double f2, double ratio,
                               double* d_coarse, double* d_fine)
{
    const double d0 = fabs(f0 - f1);
    const double d1 = fabs(f1 - f2);
    if (d_coarse) *d_coarse = d0;
    if (d_fine)   *d_fine   = d1;
    if (!(d0 > 0.0) || !(d1 > 0.0) || !(ratio > 1.0)) return NAN;
    return log(d0 / d1) / log(ratio);
}

/* One row of the unified Richardson table. */
static void write_richardson_rows(FILE* f, const char* sweep, const char* integrator,
                                  const char* refined_quantity, double ratio,
                                  int photoperiod, const Metrics* m, int nlev,
                                  const double* level_value)
{
    for (int k = 0; k < NMETRIC; ++k) {
        for (int i = 0; i + 2 < nlev; ++i) {
            const double f0 = metric_get(&m[i],     k);
            const double f1 = metric_get(&m[i + 1], k);
            const double f2 = metric_get(&m[i + 2], k);
            double d0 = 0.0, d1 = 0.0;
            const double po = richardson_order(f0, f1, f2, ratio, &d0, &d1);
            fprintf(f, "%s,%s,%s,%.10g,%d,%s,%.10g,%.10g,%.10g,"
                       "%.12g,%.12g,%.12g,%.6e,%.6e,%.6f\n",
                sweep, integrator, refined_quantity, ratio,
                photoperiod, METRIC_NAME[k],
                level_value[i], level_value[i + 1], level_value[i + 2],
                f0, f1, f2, d0, d1, po);
        }
    }
}

/* ------------------------------------------------------------------ main */

int main(int argc, char** argv)
{
    const char* config_path = "config/default.config";
    const char* out_dir     = "results/convergence";
    int         days_override = 0;

    for (int i = 1; i < argc; ++i) {
        if      (strcmp(argv[i], "--config") == 0 && i + 1 < argc) config_path   = argv[++i];
        else if (strcmp(argv[i], "--out")    == 0 && i + 1 < argc) out_dir       = argv[++i];
        else if (strcmp(argv[i], "--days")   == 0 && i + 1 < argc) days_override = atoi(argv[++i]);
        else {
            fprintf(stderr, "Usage: %s [--config <path>] [--out <dir>] [--days N]\n", argv[0]);
            return 1;
        }
    }

    /* config_load_file() assigns only the keys present in the file, so both
     * structs must carry their defaults first; otherwise every unlisted field
     * retains uninitialised stack contents and a different model is silently
     * integrated. */
    SimulationConfig cfg   = simulation_default_config();
    Input            input = generate_input();
    if (config_load_file(config_path, &cfg, &input) != 0) {
        fprintf(stderr, "Failed to load config: %s\n", config_path);
        return 1;
    }
    cfg.gui_enabled = 0;
    cfg.write_csv   = 0;
    if (days_override > 0) cfg.days = days_override;

    const real_t dt_prod  = cfg.dt_hours;
    const int    days     = cfg.days;
    const real_t psi_prod = input.gas_exchange.leaf_water_potential;

    real_t eps_ci_prod, eps_an_prod;
    ci_solver_get_tolerances(&eps_ci_prod, &eps_an_prod);

    ensure_dir("results");
    ensure_dir(out_dir);

    printf("=== Numerical convergence study (reviewer issue #15) ===\n\n");
    printf("Production settings read from %s:\n", config_path);
    printf("  days                = %d\n", days);
    printf("  dt_hours            = %g  (%ld steps/day)\n",
           (double)dt_prod, (long)llround(24.0 / (double)dt_prod));
    printf("  solver_type         = %s\n",
           cfg.solver_type == SOLVER_RK4 ? "RK4"
           : cfg.solver_type == SOLVER_ODE45 ? "ODE45" : "RKF78");
    printf("  solver_tolerance    = %g\n", (double)cfg.solver_tolerance);
    printf("  Ci solver mode      = %s\n",
           ci_solver_get_mode() == CI_SOLVER_ANALYTIC ? "analytic (closed form)"
                                                      : "fixed point");
    printf("  Ci eps_ci / eps_an  = %g / %g  (eps_an < 0 => eps_ci applied to both)\n",
           (double)eps_ci_prod, (double)eps_an_prod);
    printf("  leaf_water_potential= %g MPa\n\n", (double)psi_prod);

    char path[512];

    /* The unified Richardson table collects the observed orders of all three
     * refinement sweeps, each row labelled with the integrator and the quantity
     * actually refined, so that a tolerance order can never be mistaken for a
     * step-size order. */
    snprintf(path, sizeof(path), "%s/convergence_richardson.csv", out_dir);
    FILE* frich = fopen(path, "w");
    if (!frich) { fprintf(stderr, "cannot write %s\n", path); return 1; }
    fprintf(frich, "sweep,integrator,refined_quantity,refinement_ratio,"
                   "photoperiod_h,metric,level_coarse,level_mid,level_fine,"
                   "f_coarse,f_mid,f_fine,diff_coarse,diff_fine,observed_order\n");

    /* ------------------------- (a) ODE45 tolerance refinement at production dt */

    /* THE PRODUCTION-RELEVANT CONVERGENCE TEST.  dt is pinned at the production
     * 0.1 h and solver_tolerance is swept, because under the adaptive
     * integrator the tolerance is what controls the discretisation error. */

    static Metrics stolm[NPP][NSTOL];

    ci_solver_set_mode(CI_SOLVER_ANALYTIC);

    printf("-- (a) PRODUCTION ODE45 solver-tolerance refinement"
           " (dt fixed at %g h, analytic Ci) --\n", (double)dt_prod);
    printf("     dt is h_init + the observer sampling interval, so the tolerance"
           " is the accuracy knob.\n");
    snprintf(path, sizeof(path), "%s/convergence_solver_tol.csv", out_dir);
    FILE* f = fopen(path, "w");
    if (!f) { fprintf(stderr, "cannot write %s\n", path); return 1; }
    fprintf(f, "photoperiod_h,days,integrator,dt_hours,solver_tolerance,"
               "is_production,diverged,"
               "rgr,final_biomass,eod_starch,eod_sucrose,eod_sucrose_total,"
               "rel_dev_rgr_vs_tightest,rel_dev_biomass_vs_tightest,"
               "rel_dev_starch_vs_tightest,rel_dev_sucrose_vs_tightest,"
               "rel_dev_sucrose_total_vs_tightest,"
               "rel_dev_rgr_vs_production,rel_dev_starch_vs_production\n");

    for (int p = 0; p < NPP; ++p) {
        for (int e = 0; e < NSTOL; ++e) {
            SimulationConfig ce = cfg;
            ce.solver_type      = SOLVER_ODE45;
            ce.dt_hours         = dt_prod;
            ce.solver_tolerance = (real_t)STOL[e];
            if (run_production(&ce, &input, &SCENARIOS[p], &stolm[p][e]) != 0) {
                fprintf(stderr, "ODE45 tolerance sweep failed: %d h, tol = %g\n",
                        SCENARIOS[p].photoperiod_h, STOL[e]);
                return 1;
            }
            if (stolm[p][e].diverged)
                fprintf(stderr, "warning: ODE45 run classified DIVERGENT at "
                                "%d h, tol = %g\n",
                        SCENARIOS[p].photoperiod_h, STOL[e]);
        }
        const Metrics* tight = &stolm[p][NSTOL - 1];
        const Metrics* prod  = &stolm[p][STOL_PROD_IDX];
        for (int e = 0; e < NSTOL; ++e) {
            const Metrics* m = &stolm[p][e];
            fprintf(f, "%d,%d,ODE45,%.10g,%.3e,%d,%d,"
                       "%.12g,%.12g,%.12g,%.12g,%.12g,"
                       "%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n",
                SCENARIOS[p].photoperiod_h, days, (double)dt_prod, STOL[e],
                (e == STOL_PROD_IDX) ? 1 : 0, m->diverged,
                m->rgr, m->final_biomass, m->eod_starch, m->eod_sucrose,
                m->eod_sucrose_total,
                rel_dev(m->rgr, tight->rgr),
                rel_dev(m->final_biomass, tight->final_biomass),
                rel_dev(m->eod_starch, tight->eod_starch),
                rel_dev(m->eod_sucrose, tight->eod_sucrose),
                rel_dev(m->eod_sucrose_total, tight->eod_sucrose_total),
                rel_dev(m->rgr, prod->rgr),
                rel_dev(m->eod_starch, prod->eod_starch));
        }
        printf("  %2dh  RGR: ", SCENARIOS[p].photoperiod_h);
        for (int e = 0; e < NSTOL; ++e) printf("%.10f ", stolm[p][e].rgr);
        printf("\n       rel dev vs tightest (RGR): ");
        for (int e = 0; e < NSTOL; ++e)
            printf("%.3e ", rel_dev(stolm[p][e].rgr, stolm[p][NSTOL - 1].rgr));
        printf("\n");

        write_richardson_rows(frich, "solver_tolerance", "ODE45",
                              "solver_tolerance", STOL_RATIO,
                              SCENARIOS[p].photoperiod_h,
                              stolm[p], NSTOL, STOL);
    }
    fclose(f);
    printf("  -> %s/convergence_solver_tol.csv\n\n", out_dir);

    /* ---------------- (b) ODE45 dt refinement at production tolerance 1e-8 */

    /* Same dt grid as sweep (c) but through the production integrator, where dt
     * is h_init plus the sampling interval rather than the integration step.
     * Whatever moves here is an output-grid effect. */

    double DTP[NDTP];
    DTP[0] = (double)dt_prod * 4.0;
    DTP[1] = (double)dt_prod * 2.0;
    DTP[2] = (double)dt_prod;
    for (int i = 3; i < NDTP; ++i) DTP[i] = DTP[i - 1] * 0.5;
    const int DTP_PROD_IDX = 2;

    static Metrics dtpm[NPP][NDTP];

    printf("-- (b) PRODUCTION ODE45 dt refinement (solver_tolerance fixed at %g) --\n",
           (double)cfg.solver_tolerance);
    printf("     dt here is h_init + sampling interval, NOT the integration step.\n");
    snprintf(path, sizeof(path), "%s/convergence_dt_ode45.csv", out_dir);
    f = fopen(path, "w");
    if (!f) { fprintf(stderr, "cannot write %s\n", path); return 1; }
    fprintf(f, "photoperiod_h,days,integrator,dt_hours,samples_per_day,"
               "solver_tolerance,is_production,diverged,dt_role,"
               "rgr,final_biomass,eod_starch,eod_sucrose,eod_sucrose_total,"
               "rel_dev_rgr_vs_finest,rel_dev_biomass_vs_finest,"
               "rel_dev_starch_vs_finest,rel_dev_sucrose_vs_finest,"
               "rel_dev_sucrose_total_vs_finest,"
               "rel_dev_rgr_vs_production,rel_dev_starch_vs_production\n");

    for (int p = 0; p < NPP; ++p) {
        for (int i = 0; i < NDTP; ++i) {
            SimulationConfig cd = cfg;
            cd.solver_type      = SOLVER_ODE45;
            cd.dt_hours         = (real_t)DTP[i];
            cd.solver_tolerance = cfg.solver_tolerance;
            if (run_production(&cd, &input, &SCENARIOS[p], &dtpm[p][i]) != 0) {
                fprintf(stderr, "ODE45 dt sweep failed: %d h, dt = %g\n",
                        SCENARIOS[p].photoperiod_h, DTP[i]);
                return 1;
            }
            if (dtpm[p][i].diverged)
                fprintf(stderr, "warning: ODE45 run classified DIVERGENT at "
                                "%d h, dt = %g\n",
                        SCENARIOS[p].photoperiod_h, DTP[i]);
        }
        const Metrics* fine = &dtpm[p][NDTP - 1];
        const Metrics* prod = &dtpm[p][DTP_PROD_IDX];
        for (int i = 0; i < NDTP; ++i) {
            const Metrics* m = &dtpm[p][i];
            fprintf(f, "%d,%d,ODE45,%.10g,%ld,%.3e,%d,%d,h_init_and_sampling,"
                       "%.12g,%.12g,%.12g,%.12g,%.12g,"
                       "%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n",
                SCENARIOS[p].photoperiod_h, days, DTP[i], m->steps_per_day,
                (double)cfg.solver_tolerance,
                (i == DTP_PROD_IDX) ? 1 : 0, m->diverged,
                m->rgr, m->final_biomass, m->eod_starch, m->eod_sucrose,
                m->eod_sucrose_total,
                rel_dev(m->rgr, fine->rgr),
                rel_dev(m->final_biomass, fine->final_biomass),
                rel_dev(m->eod_starch, fine->eod_starch),
                rel_dev(m->eod_sucrose, fine->eod_sucrose),
                rel_dev(m->eod_sucrose_total, fine->eod_sucrose_total),
                rel_dev(m->rgr, prod->rgr),
                rel_dev(m->eod_starch, prod->eod_starch));
        }
        printf("  %2dh  RGR: ", SCENARIOS[p].photoperiod_h);
        for (int i = 0; i < NDTP; ++i) printf("%.10f ", dtpm[p][i].rgr);
        printf("\n       rel dev vs finest (RGR): ");
        for (int i = 0; i < NDTP; ++i)
            printf("%.3e ", rel_dev(dtpm[p][i].rgr, dtpm[p][NDTP - 1].rgr));
        printf("\n");

        write_richardson_rows(frich, "ode45_dt", "ODE45",
                              "dt_hours_h_init_and_sampling", 2.0,
                              SCENARIOS[p].photoperiod_h,
                              dtpm[p], NDTP, DTP);
    }
    fclose(f);
    printf("  -> %s/convergence_dt_ode45.csv\n\n", out_dir);

    /* ------------------------------------------- error-source separation table */

    /* The single table the manuscript needs: at the production point, how much
     * of the reported number is governed by the tolerance and how much by the
     * output grid.  Both columns are spans over the swept range with the other
     * knob held at production, so they are directly comparable. */

    printf("-- error-source separation at the production point --\n");
    snprintf(path, sizeof(path), "%s/convergence_error_separation.csv", out_dir);
    f = fopen(path, "w");
    if (!f) { fprintf(stderr, "cannot write %s\n", path); return 1; }
    fprintf(f, "photoperiod_h,days,integrator,metric,production_value,"
               "tol_span_rel_dev,tol_span_from,tol_span_to,dt_fixed_h,"
               "dt_span_rel_dev,dt_span_from,dt_span_to,solver_tolerance_fixed,"
               "dominant_error_source,tol_over_dt_ratio\n");

    for (int p = 0; p < NPP; ++p) {
        for (int k = 0; k < NMETRIC; ++k) {
            /* Span of the tolerance sweep: coarsest level against tightest. */
            const double tol_span = rel_dev(metric_get(&stolm[p][0], k),
                                            metric_get(&stolm[p][NSTOL - 1], k));
            /* Span of the dt sweep: coarsest level against finest. */
            const double dt_span  = rel_dev(metric_get(&dtpm[p][0], k),
                                            metric_get(&dtpm[p][NDTP - 1], k));
            const double ratio = (dt_span > 0.0) ? tol_span / dt_span : INFINITY;
            fprintf(f, "%d,%d,ODE45,%s,%.12g,%.6e,%.3e,%.3e,%.10g,"
                       "%.6e,%.10g,%.10g,%.3e,%s,%.6g\n",
                SCENARIOS[p].photoperiod_h, days, METRIC_NAME[k],
                metric_get(&stolm[p][STOL_PROD_IDX], k),
                tol_span, STOL[0], STOL[NSTOL - 1], (double)dt_prod,
                dt_span, DTP[0], DTP[NDTP - 1], (double)cfg.solver_tolerance,
                (tol_span > dt_span) ? "solver_tolerance" : "dt_output_grid",
                ratio);
            if (k == 0)
                printf("  %2dh  RGR: tolerance span %.3e   dt span %.3e"
                       "   ratio %.4g  -> %s dominates\n",
                    SCENARIOS[p].photoperiod_h, tol_span, dt_span, ratio,
                    (tol_span > dt_span) ? "TOLERANCE" : "dt");
        }
    }
    fclose(f);
    printf("  -> %s/convergence_error_separation.csv\n\n", out_dir);

    /* ---------------------------------------------------------------- (c) dt */

    /* Seven levels, production value bracketed from above and below.  All
     * divide 24 and every photoperiod in {4,6,8,12} exactly. */
    real_t DT[NDT];
    DT[0] = dt_prod * REAL(4.0);
    DT[1] = dt_prod * REAL(2.0);
    DT[2] = dt_prod;
    for (int i = 3; i < NDT; ++i) DT[i] = DT[i - 1] * REAL(0.5);

    static Metrics dtm[NPP][NDT];

    ci_solver_set_mode(CI_SOLVER_ANALYTIC);

    printf("-- (c) FIXED-STEP RK4 step-size refinement (analytic Ci solver) --\n");
    printf("     BACKEND VERIFICATION ONLY. Production is ODE45; no number in\n"
           "     this sweep describes a published run. dt here IS the\n"
           "     integration step, which is why the observed order is\n"
           "     meaningful and the carbon balance is exact.\n");
    snprintf(path, sizeof(path), "%s/convergence_dt.csv", out_dir);
    f = fopen(path, "w");
    if (!f) { fprintf(stderr, "cannot write %s\n", path); return 1; }
    /* Columns appended only, never inserted or reordered: is_production_config
     * and integrator are new tail columns that correct the old header's
     * implication that the production run lives on this sweep.  is_production
     * is kept, with its original meaning (dt equals the production dt_hours),
     * so existing readers keep working. */
    fprintf(f, "photoperiod_h,days,dt_hours,steps_per_day,is_production,"
               "rgr,final_biomass,eod_starch,eod_sucrose,eod_sucrose_total,"
               "int_assimilation,carbon_residual_abs,carbon_residual_rel,"
               "rel_dev_rgr_vs_finest,rel_dev_biomass_vs_finest,"
               "rel_dev_starch_vs_finest,rel_dev_sucrose_vs_finest,"
               "rel_dev_sucrose_total_vs_finest,"
               "integrator,is_production_config,dt_role\n");

    for (int p = 0; p < NPP; ++p) {
        for (int i = 0; i < NDT; ++i) {
            if (run_rk4(&input, &SCENARIOS[p], days, DT[i], &dtm[p][i]) != 0) {
                fprintf(stderr, "dt sweep failed: %d h, dt = %g\n",
                        SCENARIOS[p].photoperiod_h, (double)DT[i]);
                return 1;
            }
        }
        const Metrics* fine = &dtm[p][NDT - 1];
        for (int i = 0; i < NDT; ++i) {
            const Metrics* m = &dtm[p][i];
            fprintf(f, "%d,%d,%.10g,%ld,%d,%.12g,%.12g,%.12g,%.12g,%.12g,"
                       "%.12g,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,"
                       "RK4_fixed_step,0,integration_step\n",
                SCENARIOS[p].photoperiod_h, days, (double)DT[i], m->steps_per_day,
                (fabs((double)DT[i] - (double)dt_prod) < 1e-15) ? 1 : 0,
                m->rgr, m->final_biomass, m->eod_starch, m->eod_sucrose,
                m->eod_sucrose_total, m->int_assimilation,
                m->carbon_residual_abs, m->carbon_residual_rel,
                rel_dev(m->rgr, fine->rgr),
                rel_dev(m->final_biomass, fine->final_biomass),
                rel_dev(m->eod_starch, fine->eod_starch),
                rel_dev(m->eod_sucrose, fine->eod_sucrose),
                rel_dev(m->eod_sucrose_total, fine->eod_sucrose_total));
        }
        printf("  %2dh  RGR: ", SCENARIOS[p].photoperiod_h);
        for (int i = 0; i < NDT; ++i) printf("%.9f ", dtm[p][i].rgr);
        printf("\n       carbon residual (rel): ");
        for (int i = 0; i < NDT; ++i) printf("%.2e ", dtm[p][i].carbon_residual_rel);
        printf("\n");
    }
    fclose(f);
    printf("  -> %s\n\n", path);

    /* ------------------------------------------- observed order of convergence */

    printf("-- (c) observed order of convergence, Richardson on successive triples --\n");
    printf("     RK4 fixed-step ONLY. This is the only sweep in which dt is the\n"
           "     integration step, so it is the only one whose order can be\n"
           "     compared against a formal method order at all.\n");
    snprintf(path, sizeof(path), "%s/convergence_dt_order.csv", out_dir);
    f = fopen(path, "w");
    if (!f) { fprintf(stderr, "cannot write %s\n", path); return 1; }
    fprintf(f, "photoperiod_h,metric,dt_coarse,dt_mid,dt_fine,"
               "f_coarse,f_mid,f_fine,diff_coarse,diff_fine,observed_order,"
               "integrator,refined_quantity,is_production_config\n");

    {   /* Same triples, into the unified table, labelled by integrator. */
        double DTD[NDT];
        for (int i = 0; i < NDT; ++i) DTD[i] = (double)DT[i];
        for (int p = 0; p < NPP; ++p)
            write_richardson_rows(frich, "rk4_dt", "RK4_fixed_step",
                                  "dt_hours_integration_step", 2.0,
                                  SCENARIOS[p].photoperiod_h,
                                  dtm[p], NDT, DTD);
    }

    for (int p = 0; p < NPP; ++p) {
        printf("  %2dh  ", SCENARIOS[p].photoperiod_h);
        for (int k = 0; k < NMETRIC; ++k) {
            double last_p = NAN;
            for (int i = 0; i + 2 < NDT; ++i) {
                const double f0 = metric_get(&dtm[p][i],     k);
                const double f1 = metric_get(&dtm[p][i + 1], k);
                const double f2 = metric_get(&dtm[p][i + 2], k);
                const double d0 = fabs(f0 - f1);
                const double d1 = fabs(f1 - f2);
                const double po = (d1 > 0.0 && d0 > 0.0) ? log2(d0 / d1) : NAN;
                if (i + 2 == NDT - 1) last_p = po;
                fprintf(f, "%d,%s,%.10g,%.10g,%.10g,%.12g,%.12g,%.12g,%.6e,%.6e,%.6f,"
                           "RK4_fixed_step,dt_hours_integration_step,0\n",
                    SCENARIOS[p].photoperiod_h, METRIC_NAME[k],
                    (double)DT[i], (double)DT[i + 1], (double)DT[i + 2],
                    f0, f1, f2, d0, d1, po);
            }
            printf("%s p=%.2f  ", METRIC_NAME[k], last_p);
        }
        printf("\n");
    }
    fclose(f);
    printf("  -> %s\n\n", path);

    /* --------------------------------------------------------- (b) tolerance */

    /* The analytic solver has no tolerance, so this sweep necessarily runs the
     * legacy fixed-point scheme; 1e-4 is its production value. */
    const double EPS[] = { 1e-3, 1e-4, 1e-6, 1e-8, 1e-10 };
    const int    NEPS  = (int)(sizeof(EPS) / sizeof(EPS[0]));

    printf("-- (d) Ci solver tolerance sweep (fixed-point mode, psi = %g MPa) --\n",
           (double)psi_prod);
    snprintf(path, sizeof(path), "%s/convergence_tolerance.csv", out_dir);
    f = fopen(path, "w");
    if (!f) { fprintf(stderr, "cannot write %s\n", path); return 1; }
    fprintf(f, "photoperiod_h,days,dt_hours,ci_mode,eps_ci,eps_an,is_production,"
               "rgr,final_biomass,eod_starch,eod_sucrose,eod_sucrose_total,"
               "carbon_residual_rel,ci_calls,ci_nonconvergent,ci_nonconvergent_pct,"
               "ci_mean_iterations,rel_dev_rgr_vs_tightest,"
               "rel_dev_biomass_vs_tightest,rel_dev_starch_vs_tightest\n");

    ci_solver_set_mode(CI_SOLVER_FIXED_POINT);
    static Metrics tolm[NPP][8];
    for (int p = 0; p < NPP; ++p) {
        for (int e = 0; e < NEPS; ++e) {
            ci_solver_set_tolerances((real_t)EPS[e], (real_t)EPS[e]);
            if (run_rk4(&input, &SCENARIOS[p], days, dt_prod, &tolm[p][e]) != 0) {
                fprintf(stderr, "tolerance sweep failed: %d h, eps = %g\n",
                        SCENARIOS[p].photoperiod_h, EPS[e]);
                return 1;
            }
        }
        const Metrics* tight = &tolm[p][NEPS - 1];
        for (int e = 0; e < NEPS; ++e) {
            const Metrics* m = &tolm[p][e];
            fprintf(f, "%d,%d,%.10g,fixed_point,%.3e,%.3e,%d,"
                       "%.12g,%.12g,%.12g,%.12g,%.12g,%.6e,%ld,%ld,%.4f,%.4f,"
                       "%.6e,%.6e,%.6e\n",
                SCENARIOS[p].photoperiod_h, days, (double)dt_prod,
                EPS[e], EPS[e], (fabs(EPS[e] - 1e-4) < 1e-30) ? 1 : 0,
                m->rgr, m->final_biomass, m->eod_starch, m->eod_sucrose,
                m->eod_sucrose_total, m->carbon_residual_rel,
                m->ci_calls, m->ci_nonconvergent,
                (m->ci_calls > 0) ? 100.0 * (double)m->ci_nonconvergent / (double)m->ci_calls : 0.0,
                m->ci_mean_iterations,
                rel_dev(m->rgr, tight->rgr),
                rel_dev(m->final_biomass, tight->final_biomass),
                rel_dev(m->eod_starch, tight->eod_starch));
        }
        printf("  %2dh  RGR: ", SCENARIOS[p].photoperiod_h);
        for (int e = 0; e < NEPS; ++e) printf("%.10f ", tolm[p][e].rgr);
        printf("\n       mean iters: ");
        for (int e = 0; e < NEPS; ++e) printf("%.2f ", tolm[p][e].ci_mean_iterations);
        printf("\n");
    }
    fclose(f);
    ci_solver_set_tolerances(eps_ci_prod, eps_an_prod);
    printf("  -> %s\n\n", path);

    /* -------------------------------------------------------- (c) solver mode */

    /* Water potentials spanning the fixed-point map's contractive regime and
     * its failure band.
     *
     * NOTE on the threshold: math_audit.md section 4.1 locates the onset of
     * non-convergence near psi = -1.15 MPa, but that census was run at the
     * single-leaf reference conditions of tools/solver_verification.c, which
     * use Vcmax = 80 umol m^-2 s^-1. Production runs use the config value
     * Vcmax = 45, and the threshold moves with it: at Vcmax = 45 the map is
     * still contractive at -1.3 MPa (46 iterations) and only fails at -1.4.
     * The failure is therefore a property of the (psi, Vcmax) pair -- higher
     * Vcmax steepens the demand curve and destabilises the fixed point -- not
     * of psi alone. -1.4 MPa is used here so that the comparison genuinely
     * exercises the non-convergent regime at PRODUCTION Vcmax. */
    const double PSI[] = { (double)psi_prod, -1.4 };
    const int    NPSI  = 2;

    printf("-- (e) analytic vs legacy fixed-point Ci solver --\n");
    snprintf(path, sizeof(path), "%s/convergence_solver_mode.csv", out_dir);
    f = fopen(path, "w");
    if (!f) { fprintf(stderr, "cannot write %s\n", path); return 1; }
    fprintf(f, "photoperiod_h,days,dt_hours,psi_leaf_MPa,is_production_psi,ci_mode,"
               "eps_ci,rgr,final_biomass,eod_starch,eod_sucrose,eod_sucrose_total,"
               "carbon_residual_rel,ci_calls,ci_nonconvergent,ci_nonconvergent_pct,"
               "ci_mean_iterations,rel_dev_rgr_vs_analytic,"
               "rel_dev_biomass_vs_analytic,rel_dev_starch_vs_analytic,"
               "rel_dev_sucrose_total_vs_analytic\n");

    for (int q = 0; q < NPSI; ++q) {
        Input in_psi = input;
        in_psi.gas_exchange.leaf_water_potential = (real_t)PSI[q];
        printf("  psi = %+.2f MPa\n", PSI[q]);

        for (int p = 0; p < NPP; ++p) {
            Metrics ma, mi;

            ci_solver_set_mode(CI_SOLVER_ANALYTIC);
            if (run_rk4(&in_psi, &SCENARIOS[p], days, dt_prod, &ma) != 0) return 1;

            ci_solver_set_mode(CI_SOLVER_FIXED_POINT);
            ci_solver_set_tolerances(REAL(1e-4), REAL(1e-4));
            if (run_rk4(&in_psi, &SCENARIOS[p], days, dt_prod, &mi) != 0) return 1;

            const Metrics* rows[2] = { &ma, &mi };
            const char*    names[2] = { "analytic", "fixed_point" };
            for (int r = 0; r < 2; ++r) {
                const Metrics* m = rows[r];
                fprintf(f, "%d,%d,%.10g,%.4f,%d,%s,%.3e,"
                           "%.12g,%.12g,%.12g,%.12g,%.12g,%.6e,%ld,%ld,%.4f,%.4f,"
                           "%.6e,%.6e,%.6e,%.6e\n",
                    SCENARIOS[p].photoperiod_h, days, (double)dt_prod, PSI[q],
                    (fabs(PSI[q] - (double)psi_prod) < 1e-12) ? 1 : 0,
                    names[r], (r == 0) ? 0.0 : 1e-4,
                    m->rgr, m->final_biomass, m->eod_starch, m->eod_sucrose,
                    m->eod_sucrose_total, m->carbon_residual_rel,
                    m->ci_calls, m->ci_nonconvergent,
                    (m->ci_calls > 0) ? 100.0 * (double)m->ci_nonconvergent / (double)m->ci_calls : 0.0,
                    m->ci_mean_iterations,
                    rel_dev(m->rgr, ma.rgr),
                    rel_dev(m->final_biomass, ma.final_biomass),
                    rel_dev(m->eod_starch, ma.eod_starch),
                    rel_dev(m->eod_sucrose_total, ma.eod_sucrose_total));
            }
            printf("    %2dh  RGR analytic %.10f  fixed-point %.10f  rel dev %.3e"
                   "   (nonconv %.2f %%)\n",
                SCENARIOS[p].photoperiod_h, ma.rgr, mi.rgr,
                rel_dev(mi.rgr, ma.rgr),
                (mi.ci_calls > 0) ? 100.0 * (double)mi.ci_nonconvergent / (double)mi.ci_calls : 0.0);
        }
    }
    fclose(f);
    ci_solver_set_mode(CI_SOLVER_ANALYTIC);
    ci_solver_set_tolerances(eps_ci_prod, eps_an_prod);
    printf("  -> %s\n\n", path);

    /* -------------------------------------------------------- (d) integrator */

    printf("-- (f) integrator backend at production settings"
           " (dt = %g h, tol = %g) --\n",
           (double)cfg.dt_hours, (double)cfg.solver_tolerance);
    snprintf(path, sizeof(path), "%s/convergence_integrator.csv", out_dir);
    f = fopen(path, "w");
    if (!f) { fprintf(stderr, "cannot write %s\n", path); return 1; }
    fprintf(f, "photoperiod_h,days,dt_hours,backend,solver_tolerance,"
               "rgr,final_biomass,eod_starch,eod_sucrose,eod_sucrose_total,"
               "rel_dev_rgr_vs_config_backend,rel_dev_biomass_vs_config_backend,"
               "rel_dev_starch_vs_config_backend\n");

    /* RKF78 is deliberately excluded. The backend in the linked solver library
     * is not a consistent method: on dy/dt = -y it returns y(t) = exp(-0.902381 t)
     * instead of exp(-t), independently of the requested tolerance from 1e-2 to
     * 1e-6, and hangs (step size driven to zero) for tol <= 1e-10. The measured
     * decay rate is 0.902380952 = 1 - 82/840 to seven figures, i.e. the stage
     * weights sum to 1 - 2*(41/840) rather than 1: the two 41/840 coefficients
     * that Fehlberg's 8th-order formula assigns to the first and eleventh stages
     * are missing from the weight vector. A method whose weights do not sum to
     * unity has order zero and cannot reproduce even a linear ODE, so including
     * it in a convergence comparison would measure the library defect rather
     * than this model's discretisation error. Only the library ships as a
     * binary archive here (no source), so it cannot be repaired from this repo.
     * Reported in docs/revision/convergence_study.md; RK4 and ODE45 both pass
     * the same linear test to their design order. */
    const SolverType BK[2]      = { SOLVER_ODE45, SOLVER_RK4 };
    const char*      BK_NAME[2] = { "ODE45", "RK4" };
    const int        NBK        = 2;

    for (int p = 0; p < NPP; ++p) {
        Metrics mref;
        SimulationConfig cref = cfg;
        if (run_production(&cref, &input, &SCENARIOS[p], &mref) != 0) return 1;

        for (int b = 0; b < NBK; ++b) {
            SimulationConfig cb = cfg;
            cb.solver_type = BK[b];
            Metrics m;
            if (run_production(&cb, &input, &SCENARIOS[p], &m) != 0) return 1;
            fprintf(f, "%d,%d,%.10g,%s,%.3e,%.12g,%.12g,%.12g,%.12g,%.12g,"
                       "%.6e,%.6e,%.6e\n",
                SCENARIOS[p].photoperiod_h, days, (double)cfg.dt_hours,
                BK_NAME[b], (double)cfg.solver_tolerance,
                m.rgr, m.final_biomass, m.eod_starch, m.eod_sucrose,
                m.eod_sucrose_total,
                rel_dev(m.rgr, mref.rgr),
                rel_dev(m.final_biomass, mref.final_biomass),
                rel_dev(m.eod_starch, mref.eod_starch));
            printf("    %2dh  %-6s RGR %.10f  rel dev vs config backend %.3e\n",
                SCENARIOS[p].photoperiod_h, BK_NAME[b], m.rgr,
                rel_dev(m.rgr, mref.rgr));
        }
    }
    fclose(f);
    printf("  -> %s\n\n", path);

    /* ------------------------------------------- (e) backend sanity, exact ODE */

    /* Sweep (f) compares backends against each other, which can only detect
     * disagreement, not error. This sweep integrates dy/dt = -y, whose exact
     * solution y(t) = exp(-t) is known, so each backend's absolute error is
     * measurable. It is what identified the RKF78 defect documented above, and
     * it is the check any future change to the solver library must pass. */

    printf("-- (g) integrator verification against dy/dt = -y, y(0) = 1 --\n");
    snprintf(path, sizeof(path), "%s/convergence_integrator_exact.csv", out_dir);
    f = fopen(path, "w");
    if (!f) { fprintf(stderr, "cannot write %s\n", path); return 1; }
    fprintf(f, "backend,tol,t_end,y_numeric,y_exact,abs_error,rel_error,"
               "implied_decay_rate,consistent\n");

    const SolverType VBK[3]      = { SOLVER_RK4, SOLVER_ODE45, SOLVER_RKF78 };
    const char*      VBK_NAME[3] = { "RK4", "ODE45", "RKF78" };
    /* RKF78 hangs for tol <= 1e-10 (step size collapses), so the verification
     * grid stops at 1e-6 -- enough to show the error is tolerance-independent. */
    const double VTOL[3] = { 1e-4, 1e-6, 1e-8 };
    const double VTEND   = 1.0;

    for (int b = 0; b < 3; ++b) {
        for (int q = 0; q < 3; ++q) {
            /* Skip the RKF78/1e-8 cell: non-terminating, not merely slow. */
            if (VBK[b] == SOLVER_RKF78 && VTOL[q] < 1e-6) {
                fprintf(f, "%s,%.1e,%.3f,,,,,,hangs\n",
                    VBK_NAME[b], VTOL[q], VTEND);
                printf("    %-6s tol %.0e  DID NOT TERMINATE (step size -> 0)\n",
                    VBK_NAME[b], VTOL[q]);
                continue;
            }
            real_t y = REAL(1.0);
            SolverRunConfig rc = {
                .type    = VBK[b],
                .t0      = REAL(0.0),
                .t_end   = (real_t)VTEND,
                .h_init  = REAL(0.1),
                .tol     = (real_t)VTOL[q],
                .plugins = NULL
            };
            solver_vector_run(&rc, &y, 1, linear_decay_rhs, NULL);

            const double yex  = exp(-VTEND);
            const double aerr = fabs((double)y - yex);
            const double rate = -log((double)y) / VTEND;
            const int    ok   = (aerr < 1e-3) ? 1 : 0;
            fprintf(f, "%s,%.1e,%.3f,%.16g,%.16g,%.6e,%.6e,%.12f,%s\n",
                VBK_NAME[b], VTOL[q], VTEND, (double)y, yex,
                aerr, aerr / yex, rate, ok ? "yes" : "NO");
            printf("    %-6s tol %.0e  y = %.12f  exact %.12f  abs err %.3e"
                   "  implied rate %.9f%s\n",
                VBK_NAME[b], VTOL[q], (double)y, yex, aerr, rate,
                ok ? "" : "   <-- INCONSISTENT");
        }
    }
    fclose(f);
    printf("  -> %s\n\n", path);

    fclose(frich);
    printf("  -> %s/convergence_richardson.csv\n\n", out_dir);

    printf("Done. Nine CSVs written to %s/\n", out_dir);
    printf("PRODUCTION-RELEVANT: convergence_solver_tol.csv (ODE45, tolerance).\n");
    printf("NOT production:      convergence_dt.csv (fixed-step RK4 verification).\n");
    return 0;
}

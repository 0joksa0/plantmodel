/*
 * carbon_audit.c -- term-by-term carbon budget of the model against the
 * published budget of Sulpice et al. (2014), supplement mmc3.xlsx sheet
 * "Growth components".
 *
 * Relation to tools/carbon_balance.c
 * ----------------------------------
 * carbon_balance.c answers "does the carbohydrate subsystem conserve carbon
 * numerically?" using a fixed-step augmented RK4.  That is a closure check and
 * the integrator does not have to match production, because a closure residual
 * is a statement about arithmetic.
 *
 * This tool answers a different question: "how large is each carbon flux, and
 * how does it compare with the measured budget?"  A flux magnitude is only
 * meaningful on the production trajectory, so the state here is advanced by
 * the same ODE45 path production uses (model_run + SOLVER_ODE45 + the
 * production tolerance), not by fixed-step RK4.
 *
 * Method and its one deviation from production
 * --------------------------------------------
 * Fluxes are accumulated as extra ODE components appended to the 11 plant
 * states, so the quadrature uses exactly the stage values the state
 * integration uses.  Appending components to the state does change what the
 * adaptive controller sees: the accumulators enter the error norm.  To keep
 * that deviation confined and measurable, the run is split:
 *
 *   phase 1   t in [0, (D-1)*24]   plain X_DIM production model, bit-identical
 *                                  to simulate_days() over that interval
 *   phase 2   t in [(D-1)*24, D*24] augmented X_DIM+F_DIM model
 *
 * Only the final day is integrated with the augmented system, and the tool
 * additionally runs the plain model over the same final day from the same
 * phase-1 end state and prints the state-wise relative difference, so the
 * size of the deviation is reported rather than assumed.  --method rk4
 * selects the fixed-step augmented RK4 of carbon_balance.c as a second,
 * independent quadrature.
 *
 * Units
 * -----
 * The model's carbon currency is umol hexose gFW^-1 h^-1 (see
 * docs/revision/sla_conversion_resolved.md).  Sulpice's "Growth components"
 * sheet reports umol glucose equivalents gFW^-1 h^-1 for rates and umol C
 * gFW^-1 for daily totals; 1 umol hexose = 6 umol C.  Both are emitted.
 *
 * Everything here is opt-in: no default behaviour of any existing target is
 * touched, and src/ is not modified.
 *
 * Output
 * ------
 *   <out>/carbon_audit_terms.csv   per-photoperiod final-24h budget, both units
 *   <out>/carbon_audit_diel_<PP>h.csv  final-24h diel trace, states + fluxes
 *   <out>/carbon_audit_saturation.csv  saturation inventory per photoperiod
 *   <out>/carbon_audit_parity.csv      ODE45 augmented vs plain final states
 */

#define _POSIX_C_SOURCE 200809L

#include "config/config.h"
#include "model/input.h"
#include "model/model.h"
#include "model/state.h"
#include "simulation/dynamics.h"
#include "simulation/observers.h"
#include "simulation/simulation.h"

/* framework/solver_plugin.h uses max_align_t, which gnu99 does not define;
 * src/simulation/simulation.c applies the same shim. */
#ifndef max_align_t
#define max_align_t long double
#define CARBON_AUDIT_MAX_ALIGN_T_SHIM
#endif
#include <framework/default_runtime.h>
#include <framework/model_buffers.h>
#include <framework/model_runner.h>
#include <framework/solver_plugin.h>
#ifdef CARBON_AUDIT_MAX_ALIGN_T_SHIM
#undef max_align_t
#undef CARBON_AUDIT_MAX_ALIGN_T_SHIM
#endif
#include <observers.h>
#include <solver.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

/* ------------------------------------------------------------------ fluxes */

typedef enum {
    F_ASSIM = 0,   /* A_n as the model uses it (net of leaf Rd during light) */
    F_RD_LEAF,     /* Rd carried in the same converted units, light hours only */
    F_TO_STARCH,   /* k_s A_n                                                 */
    F_TO_SUCROSE,  /* (1-k_s) A_n                                             */
    F_DEGRADE,     /* D, starch -> sucrose                                    */
    F_UPTAKE,      /* C_up, nutrient uptake cost                              */
    F_TRANSPORT,   /* C_tr, phloem transport cost                             */
    F_MAINT,       /* r_m * s, maintenance respiration                        */
    F_LOAD,        /* phi_load * E_night, export to the growth sink           */
    F_VAC_IN,      /* k_in  * s                                               */
    F_VAC_OUT,     /* k_out * s_v                                             */
    F_GROWTH_L,    /* dB_l/dt   [gFW/h]                                       */
    F_GROWTH_R,    /* dB_r/dt   [gFW/h]                                       */
    F_LIGHT,       /* light indicator; integrates to photoperiod hours        */
    F_ENIGHT,      /* E_night                                                 */
    F_SMAX,        /* s/(s+max_sucrose)                                       */
    F_SMIN,        /* min_sucrose/(min_sucrose+s)                             */
    F_ECLAMP,      /* 1 while s > max_sucrose, i.e. E_night clamp active      */
    F_SUC,         /* cytosolic sucrose, for a time mean                      */
    F_STARCH,      /* starch, for a time mean                                 */
    F_LIMPHOTO,    /* limitation_of_photosynthetic_rate (unused in FJ mode)   */
    F_UPTAKE_CONS, /* uptake cost recomputed from the CURRENT uptake rates    */
    F_NUP,         /* nitrogen_uptake                                         */
    F_PUP,         /* phosphorus_uptake                                       */
    F_DIM
} FluxIndex;

static const char* FLUX_NAME[F_DIM] = {
    "assimilation_net", "leaf_Rd_light", "to_starch", "to_sucrose",
    "starch_degradation", "uptake_cost", "transport_cost", "maintenance_resp",
    "phloem_load", "vacuole_in", "vacuole_out", "growth_leaf_gFW",
    "growth_root_gFW", "light_hours", "E_night", "s_max_fract", "s_min_fract",
    "E_night_clamped", "sucrose_mean", "starch_mean", "starch_feedback",
    "uptake_cost_consistent", "nitrogen_uptake", "phosphorus_uptake"
};

#define AUG_DIM (X_DIM + F_DIM)

/* ------------------------------------------------------------------ rhs */

static void fill_fluxes(real_t t, const real_t* x, real_t* f, PlantCtx* ctx)
{
    Input s = *(ctx->base);
    simulation_unpack_state(&s, x);
    s.growth.total_biomass = s.growth.leaf_biomass + s.growth.root_biomass;
    simulation_compute_algebraic(&s, fmod(t, REAL(24.0)), ctx->starch_night_start);

    const real_t A_n = s.photo.photosynthesis;
    const real_t k_s = s.carbohydrates.starch_partition_coeff;
    const real_t D   = s.carbohydrates.starch_degradation_rate;
    const real_t suc = s.carbohydrates.sucrose;
    const real_t vac = s.carbohydrates.sucrose_vacuole;
    const real_t smx = s.carbohydrates.max_sucrose;
    const real_t smn = s.carbohydrates.min_sucrose;

    /* Same conversion convert_to_photosynthesis_per_gram_per_hour() applies,
     * so Rd is reported in the model's carbon currency rather than in
     * umol CO2 m^-2 s^-1. */
    const real_t kappa = REAL(1.10) * (RMAX(s.gas_exchange.specific_leaf_area, REAL(1e-6)) / REAL(0.025));

    f[F_ASSIM]      = A_n;
    f[F_RD_LEAF]    = (s.photo.light > REAL(0.0)) ? s.gas_exchange.respiration_rate * kappa : REAL(0.0);
    f[F_TO_STARCH]  = k_s * A_n;
    f[F_TO_SUCROSE] = (REAL(1.0) - k_s) * A_n;
    f[F_DEGRADE]    = D;
    f[F_UPTAKE]     = s.carbohydrates.uptake_cost;
    f[F_TRANSPORT]  = s.carbohydrates.transport_cost;
    f[F_MAINT]      = s.carbohydrates.respiration_frequency * suc;
    f[F_LOAD]       = s.carbohydrates.sucrose_loading_frequency * s.carbohydrates.night_efficiency_starch;
    f[F_VAC_IN]     = s.carbohydrates.k_vacuole_in  * suc;
    f[F_VAC_OUT]    = s.carbohydrates.k_vacuole_out * vac;
    f[F_LIGHT]      = s.photo.light;
    f[F_ENIGHT]     = s.carbohydrates.night_efficiency_starch;
    f[F_SMAX]       = suc / (suc + smx);
    f[F_SMIN]       = smn / (smn + suc);
    f[F_ECLAMP]     = (suc > smx) ? REAL(1.0) : REAL(0.0);
    f[F_SUC]        = suc;
    f[F_STARCH]     = s.carbohydrates.starch;
    f[F_LIMPHOTO]   = s.photo.limitation_of_photosynthetic_rate;

    /* plant_compute_outputs() evaluates uptake_cost at algebraic.c:86 from
     * nutrients.nitrogen_uptake / phosphorus_uptake, which are themselves only
     * written 41 lines later (algebraic.c:127, :131).  It therefore always
     * reads the config-supplied 0.0.  This diagnostic recomputes the same
     * expression from the uptake rates the SAME call produced, so the audit can
     * report the term's intended magnitude alongside the zero the model uses.
     * Read-only; src/ is unchanged. */
    /* Growth is evaluated here rather than copied from dx so that the diel
     * observer, which calls fill_fluxes() alone, records it too. */
    f[F_GROWTH_L] = leaf_growth(s.growth.lambda_sb, s.growth.sucrose_root_allocation,
                                s.carbohydrates.sucrose_loading_frequency,
                                s.carbohydrates.night_efficiency_starch,
                                s.growth.leaf_biomass, s.growth.leaf_deathrate,
                                s.growth.leaf_competitive_rate);
    f[F_GROWTH_R] = root_growth(s.growth.lambda_sb, s.growth.sucrose_root_allocation,
                                s.carbohydrates.sucrose_loading_frequency,
                                s.carbohydrates.night_efficiency_starch,
                                s.growth.leaf_biomass, s.growth.root_biomass,
                                s.growth.root_deathrate, s.growth.root_competitive_rate);

    f[F_NUP]        = s.nutrients.nitrogen_uptake;
    f[F_PUP]        = s.nutrients.phosphorus_uptake;
    f[F_UPTAKE_CONS] = uptake_cost(
        s.nutrients.nitrogen_uptake_sucrose_consumption, s.nutrients.nitrogen_uptake,
        s.nutrients.phosphorus_uptake_sucrose_consumption, s.nutrients.phosphorus_uptake);
}

static void augmented_rhs_raw(real_t t, const real_t* x, real_t* dx, PlantCtx* ctx)
{
    simulation_plant_rhs(t, x, dx, ctx);
    fill_fluxes(t, x, dx + X_DIM, ctx);
}

static void augmented_model_rhs(
    real_t t, const real_t* state, const real_t* parameters,
    const real_t* inputs, real_t* dstate, void* model_ctx)
{
    (void)parameters; (void)inputs;
    augmented_rhs_raw(t, state, dstate, (PlantCtx*)model_ctx);
}

/* model_buffers_init() rejects a descriptor without a state-field table, so
 * the 11 plant states are named and the F_DIM accumulators appended. */
#define PS(i, n) { .name = n, .kind = MODEL_FIELD_STATE, .index = (i), .unit = NULL, .description = NULL }
#define FS(i, n) { .name = "acc_" n, .kind = MODEL_FIELD_STATE, .index = X_DIM + (i), .unit = NULL, .description = NULL }

static const ModelFieldDescriptor AUG_FIELDS[AUG_DIM] = {
    PS(X_STARCH_PART, "starch_partition_coeff"), PS(X_STARCH, "starch"),
    PS(X_SUCROSE, "sucrose"), PS(X_N_AFF, "nitrogen_affinity"),
    PS(X_P_AFF, "phosphorus_affinity"), PS(X_N, "nitrogen"),
    PS(X_P, "phosphorus"), PS(X_SUCROSE_ROOT_ALLOC, "sucrose_root_allocation"),
    PS(X_LEAF_BIOMASS, "leaf_biomass"), PS(X_ROOT_BIOMASS, "root_biomass"),
    PS(X_SUCROSE_VAC, "sucrose_vacuole"),
    FS(F_ASSIM, "assimilation_net"), FS(F_RD_LEAF, "leaf_Rd_light"),
    FS(F_TO_STARCH, "to_starch"), FS(F_TO_SUCROSE, "to_sucrose"),
    FS(F_DEGRADE, "starch_degradation"), FS(F_UPTAKE, "uptake_cost"),
    FS(F_TRANSPORT, "transport_cost"), FS(F_MAINT, "maintenance_resp"),
    FS(F_LOAD, "phloem_load"), FS(F_VAC_IN, "vacuole_in"),
    FS(F_VAC_OUT, "vacuole_out"), FS(F_GROWTH_L, "growth_leaf_gFW"),
    FS(F_GROWTH_R, "growth_root_gFW"), FS(F_LIGHT, "light_hours"),
    FS(F_ENIGHT, "E_night"), FS(F_SMAX, "s_max_fract"),
    FS(F_SMIN, "s_min_fract"), FS(F_ECLAMP, "E_night_clamped"),
    FS(F_SUC, "sucrose_mean"), FS(F_STARCH, "starch_mean"),
    FS(F_LIMPHOTO, "starch_feedback"), FS(F_UPTAKE_CONS, "uptake_cost_consistent"),
    FS(F_NUP, "nitrogen_uptake"), FS(F_PUP, "phosphorus_uptake")
};
#undef PS
#undef FS

static const ModelDescriptor AUG_DESCRIPTOR = {
    .name = "plantmodel_carbon_audit",
    .state_count = AUG_DIM,
    .parameter_count = 0, .input_count = 0, .output_count = 0,
    .state_fields = AUG_FIELDS, .parameter_fields = NULL,
    .input_fields = NULL, .output_fields = NULL
};

static const ModelInterface AUG_INTERFACE = {
    .descriptor = &AUG_DESCRIPTOR,
    .callbacks = { .rhs = augmented_model_rhs, .outputs = NULL },
    .model_ctx = NULL
};

/* ------------------------------------------------------------------ observers */

typedef struct {
    FILE*    f;
    PlantCtx* ctx;
    real_t   next_t;
    real_t   dt;
    real_t   t_origin;
} DielObsCtx;

static void diel_observer(real_t t, const real_t* x, size_t n, void* vctx)
{
    (void)n;
    DielObsCtx* o = (DielObsCtx*)vctx;
    if (!o->f || t + REAL(1e-9) < o->next_t)
        return;
    o->next_t += o->dt;

    real_t f[F_DIM];
    fill_fluxes(t, x, f, o->ctx);

    fprintf(o->f, "%.6f,%.6f", (double)t, (double)(t - o->t_origin));
    for (int i = 0; i < X_DIM; ++i) fprintf(o->f, ",%.10g", (double)x[i]);
    for (int i = 0; i < F_DIM; ++i) fprintf(o->f, ",%.10g", (double)f[i]);
    fprintf(o->f, "\n");
}

/* ------------------------------------------------------------------ fixed-step RK4 (cross-check) */

static void rk4_advance(real_t t0, real_t t_end, real_t h, real_t* x,
                        PlantCtx* ctx, const Input* base, size_t dim)
{
    /* augmented_rhs_raw always writes AUG_DIM components, so the stage buffers
     * are AUG_DIM wide even when only the first `dim` are advanced (phase 1
     * updates the plant states alone and leaves the accumulators at zero). */
    real_t k1[AUG_DIM], k2[AUG_DIM], k3[AUG_DIM], k4[AUG_DIM], tp[AUG_DIM];
    memcpy(tp, x, AUG_DIM * sizeof(real_t));
    long nsteps = (long)llround((double)(t_end - t0) / (double)h);
    for (long k = 0; k < nsteps; ++k) {
        real_t t = t0 + (real_t)k * h;

        Input probe = *base;
        update_light_conditions(&probe, (real_t)fmod((double)t, 24.0));
        int is_light = (probe.photo.light != REAL(0.0));
        if (ctx->was_light && !is_light) ctx->starch_night_start = x[X_STARCH];
        ctx->was_light = is_light;

        augmented_rhs_raw(t, x, k1, ctx);
        for (size_t i = 0; i < dim; ++i) tp[i] = x[i] + REAL(0.5) * h * k1[i];
        augmented_rhs_raw(t + REAL(0.5) * h, tp, k2, ctx);
        for (size_t i = 0; i < dim; ++i) tp[i] = x[i] + REAL(0.5) * h * k2[i];
        augmented_rhs_raw(t + REAL(0.5) * h, tp, k3, ctx);
        for (size_t i = 0; i < dim; ++i) tp[i] = x[i] + h * k3[i];
        augmented_rhs_raw(t + h, tp, k4, ctx);
        for (size_t i = 0; i < dim; ++i)
            x[i] += (h / REAL(6.0)) * (k1[i] + REAL(2.0) * k2[i] + REAL(2.0) * k3[i] + k4[i]);
    }
}

/* ------------------------------------------------------------------ ODE45 driver */

/* simulate_days() installs a step budget so a collapsed adaptive step becomes
 * a classified divergence rather than a hang.  This tool drives model_run()
 * directly and must install its own, or a stalled run never returns. */
typedef struct { unsigned long steps, budget; int exhausted; } BudgetCtx;

static SolverPluginResult budget_on_attempt(void* vctx, const SolverEvent* ev)
{
    (void)ev;
    BudgetCtx* b = (BudgetCtx*)vctx;
    if (!b) return SOLVER_PLUGIN_CONTINUE;
    if (++b->steps >= b->budget) { b->exhausted = 1; return SOLVER_PLUGIN_STOP; }
    return SOLVER_PLUGIN_CONTINUE;
}

typedef struct {
    NightStartObserverCtx night_ctx;
    NanGuardCtx           nan_ctx;
    DielObsCtx            diel_ctx;
    BudgetCtx             budget_ctx;
    VectorObserver        obs[3];
    void*                 obs_ctx[3];
    ObserverChain         chain;
    SolverStats           stats;
    SolverPluginManager   plugins;
} AuditRuntime;

static void audit_runtime_init(AuditRuntime* rt, PlantCtx* ctx, FILE* diel,
                               real_t t_origin, real_t dt)
{
    memset(rt, 0, sizeof(*rt));
    rt->night_ctx = (NightStartObserverCtx){ .plant_ctx = ctx, .starch_index = X_STARCH };
    rt->diel_ctx  = (DielObsCtx){ .f = diel, .ctx = ctx, .next_t = t_origin,
                                  .dt = dt, .t_origin = t_origin };
    rt->obs[0] = nan_guard_observer;   rt->obs_ctx[0] = &rt->nan_ctx;
    rt->obs[1] = night_start_observer; rt->obs_ctx[1] = &rt->night_ctx;
    rt->obs[2] = diel ? diel_observer : NULL;
    rt->obs_ctx[2] = diel ? (void*)&rt->diel_ctx : NULL;
    rt->chain = (ObserverChain){ .obs = rt->obs, .ctx = rt->obs_ctx, .count = 3 };

    solver_stats_init(&rt->stats);
    solver_plugin_manager_init(&rt->plugins);
    SolverDefaultRuntimeConfig cfg = {
        .observer = observer_chain_cb, .observer_ctx = &rt->chain,
        .stats = &rt->stats, .logger = NULL
    };
    if (solver_default_runtime_init(&rt->plugins, &cfg) != 0)
        rt->plugins.count = 0;

    rt->budget_ctx = (BudgetCtx){ .steps = 0, .budget = 2000000u, .exhausted = 0 };
    SolverPlugin budget = { .name = "carbon_audit_step_budget",
                            .plugin_ctx = &rt->budget_ctx,
                            .on_run_start = NULL, .on_step_attempt = budget_on_attempt,
                            .on_step_accepted = NULL, .on_step_rejected = NULL,
                            .on_run_finish = NULL };
    solver_plugin_manager_add(&rt->plugins, &budget);
}

/* Advances `x` (dim components) from t0 to t_end with ODE45 at the production
 * tolerance, using `iface` (plain or augmented). */
static int ode45_advance(const ModelInterface* iface, PlantCtx* ctx,
                         real_t* x, size_t dim, real_t t0, real_t t_end,
                         real_t h_init, real_t tol, FILE* diel, real_t dt_sample,
                         int* out_nonfinite)
{
    AuditRuntime rt;
    audit_runtime_init(&rt, ctx, diel, t0, dt_sample);

    ModelInterface m = *iface;
    m.model_ctx = ctx;
    ModelBuffers buf = {0};
    if (model_buffers_init(&buf, m.descriptor) != MODEL_STATUS_OK) return 1;
    memcpy(buf.state, x, dim * sizeof(real_t));

    SolverRunConfig run = { .type = SOLVER_ODE45, .t0 = t0, .t_end = t_end,
                            .h_init = h_init, .tol = tol, .plugins = &rt.plugins };
    ModelStatus st = model_run(&m, &buf, &run);
    if (st != MODEL_STATUS_OK) { model_buffers_free(&buf); return 1; }
    memcpy(x, buf.state, dim * sizeof(real_t));
    model_buffers_free(&buf);
    if (out_nonfinite) *out_nonfinite = rt.nan_ctx.tripped || rt.budget_ctx.exhausted;
    if (rt.budget_ctx.exhausted)
        fprintf(stderr, "warning: step budget exhausted on [%g, %g]; run is divergent\n",
                (double)t0, (double)t_end);
    return 0;
}

/* ------------------------------------------------------------------ per-run record */

typedef struct {
    int    photoperiod_h;
    double flux[F_DIM];
    double a0, a1, s0, s1, v0, v1, bl0, bl1, br0, br1;
    double kappa;
    double parity_max_rel;
    int    nonfinite;
    double res_starch, res_sucrose, res_total;
} AuditRow;

static void ensure_dir(const char* p) { mkdir(p, 0775); }

/* ------------------------------------------------------------------ main */

static void usage(const char* a0)
{
    printf("Usage: %s [--config p] [--out dir] [--days N] [--kappa K] [--chamber]\n"
           "          [--sdr X] [--sdi X] [--sni X] [--method ode45|rk4] [--dt h]\n"
           "  --chamber  light_profile=1, max_PAR=160, T_day=21, T_night=19\n"
           "  --kappa    target A_n conversion factor; applied as SLA = 0.025*K/1.10\n", a0);
}

int main(int argc, char** argv)
{
    const char* config_path = "config/default.config";
    const char* out_dir     = "results/carbon_audit";
    int    days      = 30;
    double kappa     = 0.0;      /* 0 => leave SLA as configured */
    int    chamber   = 0;
    int    use_rk4   = 0;
    int    scenario  = 0;
    double dt_hours  = 0.0;
    double sdr = -1.0, sdi = -1.0, sni = -1.0;
    double max_suc = -1.0, phi_load = -1.0, r_m = -1.0, lam_g = -1.0;

    for (int i = 1; i < argc; ++i) {
        if      (!strcmp(argv[i], "--config") && i + 1 < argc) config_path = argv[++i];
        else if (!strcmp(argv[i], "--out")    && i + 1 < argc) out_dir     = argv[++i];
        else if (!strcmp(argv[i], "--days")   && i + 1 < argc) days        = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--kappa")  && i + 1 < argc) kappa       = atof(argv[++i]);
        else if (!strcmp(argv[i], "--chamber"))                chamber     = 1;
        else if (!strcmp(argv[i], "--scenario"))               scenario    = 1;
        else if (!strcmp(argv[i], "--sdr")    && i + 1 < argc) sdr         = atof(argv[++i]);
        else if (!strcmp(argv[i], "--sdi")    && i + 1 < argc) sdi         = atof(argv[++i]);
        else if (!strcmp(argv[i], "--sni")    && i + 1 < argc) sni         = atof(argv[++i]);
        else if (!strcmp(argv[i], "--dt")     && i + 1 < argc) dt_hours    = atof(argv[++i]);
        else if (!strcmp(argv[i], "--max-sucrose") && i + 1 < argc) max_suc  = atof(argv[++i]);
        else if (!strcmp(argv[i], "--phi-load")    && i + 1 < argc) phi_load = atof(argv[++i]);
        else if (!strcmp(argv[i], "--rm")          && i + 1 < argc) r_m      = atof(argv[++i]);
        else if (!strcmp(argv[i], "--lambda-g")    && i + 1 < argc) lam_g    = atof(argv[++i]);
        else if (!strcmp(argv[i], "--method") && i + 1 < argc) use_rk4     = !strcmp(argv[++i], "rk4");
        else { usage(argv[0]); return 1; }
    }

    /* config_load_file assigns only the keys present in the file. */
    SimulationConfig cfg = simulation_default_config();
    Input base           = generate_input();
    if (config_load_file(config_path, &cfg, &base) != 0) {
        fprintf(stderr, "cannot load %s\n", config_path);
        return 1;
    }
    cfg.gui_enabled = 0;
    cfg.write_csv   = 0;
    if (days > 0) cfg.days = days;
    real_t dt = (dt_hours > 0.0) ? (real_t)dt_hours : cfg.dt_hours;
    if (!(dt > REAL(0.0))) dt = REAL(0.1);

    if (chamber) {
        base.photo.light_profile = REAL(1.0);
        base.photo.max_PAR       = REAL(160.0);
        base.gas_exchange.leaf_temperature_day   = REAL(21.0);
        base.gas_exchange.leaf_temperature_night = REAL(19.0);
    }
    if (kappa > 0.0)
        base.gas_exchange.specific_leaf_area = (real_t)(0.025 * kappa / 1.10);
    if (sdr >= 0.0) base.carbohydrates.lambda_sdr = (real_t)sdr;
    if (sdi >= 0.0) base.carbohydrates.lambda_sdi = (real_t)sdi;
    if (sni >= 0.0) base.carbohydrates.lambda_sni = (real_t)sni;
    if (max_suc  > 0.0) base.carbohydrates.max_sucrose               = (real_t)max_suc;
    if (phi_load > 0.0) base.carbohydrates.sucrose_loading_frequency = (real_t)phi_load;
    if (r_m      > 0.0) base.carbohydrates.respiration_frequency     = (real_t)r_m;
    if (lam_g    > 0.0) base.carbohydrates.lambda_g                  = (real_t)lam_g;

    const double kappa_eff = 1.10 * (double)base.gas_exchange.specific_leaf_area / 0.025;

    ensure_dir("results");
    ensure_dir(out_dir);

    printf("=== carbon audit ===\n");
    printf("  config %s   days %d   dt %g   method %s\n",
           config_path, cfg.days, (double)dt, use_rk4 ? "rk4" : "ode45");
    printf("  chamber=%d  SLA=%.8g  kappa_eff=%.6f  PAR=%g  profile=%g  Td=%g Tn=%g\n",
           chamber, (double)base.gas_exchange.specific_leaf_area, kappa_eff,
           (double)base.photo.max_PAR, (double)base.photo.light_profile,
           (double)base.gas_exchange.leaf_temperature_day,
           (double)base.gas_exchange.leaf_temperature_night);
    printf("  lambda_sdr=%g sdi=%g sni=%g  max_sucrose=%g min_sucrose=%g\n\n",
           (double)base.carbohydrates.lambda_sdr, (double)base.carbohydrates.lambda_sdi,
           (double)base.carbohydrates.lambda_sni, (double)base.carbohydrates.max_sucrose,
           (double)base.carbohydrates.min_sucrose);

    const int PPS[4] = { 4, 6, 8, 12 };
    /* The per-photoperiod values tools/validation_runner.c uses (its SCENARIOS
     * table).  --scenario reproduces the published calibration; without it the
     * run uses one value per parameter for all four photoperiods. */
    const real_t SC_SNI[4] = { REAL(0.16), REAL(0.15), REAL(0.13), REAL(0.08) };
    const real_t SC_SB[4]  = { REAL(0.003843), REAL(0.004433), REAL(0.005246), REAL(0.005850) };
    const real_t SC_FB[4]  = { REAL(0.82), REAL(0.79), REAL(0.67), REAL(0.62) };

    AuditRow rows[4];
    memset(rows, 0, sizeof(rows));

    for (int p = 0; p < 4; ++p) {
        Input in = base;
        in.photo.photoperiod = (real_t)PPS[p];
        if (scenario) {
            in.carbohydrates.lambda_sni = SC_SNI[p];
            in.growth.lambda_sb         = SC_SB[p];
            in.photo.feedback_on_photosynthesis = SC_FB[p];
        }
        in.growth.lambda_sb  = lambda_sb_f(in.growth.lambda_sb,
                                           in.nutrients.nitrogen_soil_content,
                                           in.nutrients.phosphorus_soil_content);

        PlantCtx ctx = { .base = &in, .starch_night_start = in.carbohydrates.starch, .was_light = 0 };
        Input probe = in;
        update_light_conditions(&probe, REAL(0.0));
        ctx.was_light = (probe.photo.light != REAL(0.0));

        real_t x[AUG_DIM];
        memset(x, 0, sizeof(x));
        { PlantState st = plant_state_from_input(&in); plant_state_pack(x, &st); }

        const real_t t_split = (real_t)((cfg.days - 1) * 24);
        const real_t t_end   = (real_t)(cfg.days * 24);
        int nonfinite = 0;

        if (use_rk4) {
            rk4_advance(REAL(0.0), t_split, dt, x, &ctx, &in, X_DIM);
        } else {
            int nf1 = 0;
            if (ode45_advance(simulation_plant_model_interface(), &ctx, x, X_DIM,
                              REAL(0.0), t_split, dt, cfg.solver_tolerance,
                              NULL, dt, &nf1) != 0) {
                fprintf(stderr, "phase 1 failed at %d h\n", PPS[p]); return 1;
            }
            nonfinite |= nf1;
        }

        /* Snapshot for the parity check, then zero the accumulators. */
        real_t x_plain[X_DIM];
        memcpy(x_plain, x, sizeof(x_plain));
        PlantCtx ctx_plain = ctx;
        Input in_plain = in;
        ctx_plain.base = &in_plain;
        memset(x + X_DIM, 0, F_DIM * sizeof(real_t));

        char diel_path[512];
        snprintf(diel_path, sizeof(diel_path), "%s/carbon_audit_diel_%dh.csv", out_dir, PPS[p]);
        FILE* fd = fopen(diel_path, "w");
        if (fd) {
            fprintf(fd, "t_h,t_rel_h,starch_part,starch,sucrose,n_aff,p_aff,nitrogen,"
                        "phosphorus,suc_root_alloc,leaf_biomass,root_biomass,sucrose_vac");
            for (int i = 0; i < F_DIM; ++i) fprintf(fd, ",%s", FLUX_NAME[i]);
            fprintf(fd, "\n");
        }

        if (use_rk4) {
            /* Fixed-step RK4 cannot write the diel trace through the solver's
             * observer chain, so sample it here at the same cadence. */
            long nsub = (long)llround(24.0 / (double)dt);
            for (long k = 0; k < nsub; ++k) {
                real_t t = t_split + (real_t)k * dt;
                if (fd) {
                    real_t f[F_DIM];
                    fill_fluxes(t, x, f, &ctx);
                    fprintf(fd, "%.6f,%.6f", (double)t, (double)(t - t_split));
                    for (int i = 0; i < X_DIM; ++i) fprintf(fd, ",%.10g", (double)x[i]);
                    for (int i = 0; i < F_DIM; ++i) fprintf(fd, ",%.10g", (double)f[i]);
                    fprintf(fd, "\n");
                }
                rk4_advance(t, t + dt, dt, x, &ctx, &in, AUG_DIM);
            }
        } else {
            int nfa = 0;
            if (ode45_advance(&AUG_INTERFACE, &ctx, x, AUG_DIM, t_split, t_end,
                              dt, cfg.solver_tolerance, fd, dt, &nfa) != 0) {
                fprintf(stderr, "phase 2 failed at %d h\n", PPS[p]); return 1;
            }
            nonfinite |= nfa;
            /* Same final day on the plain production model, from the same
             * start state, to measure what the augmentation cost. */
            int nf2 = 0;
            ode45_advance(simulation_plant_model_interface(), &ctx_plain, x_plain, X_DIM,
                          t_split, t_end, dt, cfg.solver_tolerance, NULL, dt, &nf2);
        }
        if (fd) fclose(fd);

        AuditRow* r = &rows[p];
        r->photoperiod_h = PPS[p];
        r->kappa = kappa_eff;
        r->nonfinite = nonfinite;
        for (int i = 0; i < F_DIM; ++i) r->flux[i] = (double)x[X_DIM + i];

        double worst = 0.0;
        if (!use_rk4) {
            for (int i = 0; i < X_DIM; ++i) {
                double d = fabs((double)x[i] - (double)x_plain[i]);
                double s = fabs((double)x_plain[i]);
                double rel = (s > 1e-12) ? d / s : d;
                if (rel > worst) worst = rel;
            }
        }
        r->parity_max_rel = worst;
        r->a1 = (double)x[X_STARCH];
        r->s1 = (double)x[X_SUCROSE];
        r->v1 = (double)x[X_SUCROSE_VAC];
        r->bl1 = (double)x[X_LEAF_BIOMASS];
        r->br1 = (double)x[X_ROOT_BIOMASS];

        printf("  %2dh  A_n(24h)=%9.3f  load=%8.3f  maint=%8.3f  uptake=%7.3f  "
               "transp=%6.3f  starch=%7.3f  suc=%7.3f  parity=%.2e\n",
               PPS[p], r->flux[F_ASSIM], r->flux[F_LOAD], r->flux[F_MAINT],
               r->flux[F_UPTAKE], r->flux[F_TRANSPORT], r->a1, r->s1,
               r->parity_max_rel);
    }

    /* ---- terms CSV ---- */
    char path[512];
    snprintf(path, sizeof(path), "%s/carbon_audit_terms.csv", out_dir);
    FILE* f = fopen(path, "w");
    if (!f) { fprintf(stderr, "cannot write %s\n", path); return 1; }
    fprintf(f, "photoperiod_h,days,kappa_eff,method");
    for (int i = 0; i < F_DIM; ++i) fprintf(f, ",int24_%s", FLUX_NAME[i]);
    fprintf(f, ",starch_end,sucrose_end,vacuole_end,sucrose_total_end,"
               "leaf_biomass_end,root_biomass_end,parity_max_rel,nonfinite\n");
    for (int p = 0; p < 4; ++p) {
        AuditRow* r = &rows[p];
        fprintf(f, "%d,%d,%.10g,%s", r->photoperiod_h, cfg.days, r->kappa,
                use_rk4 ? "rk4" : "ode45");
        for (int i = 0; i < F_DIM; ++i) fprintf(f, ",%.12g", r->flux[i]);
        fprintf(f, ",%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.6e,%d\n",
                r->a1, r->s1, r->v1, r->s1 + r->v1, r->bl1, r->br1,
                r->parity_max_rel, r->nonfinite);
    }
    fclose(f);
    printf("\n  -> %s\n", path);
    return 0;
}

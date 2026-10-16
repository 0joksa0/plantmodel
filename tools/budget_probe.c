/* budget_probe -- a carbon budget that closes, measured per step.
 *
 * WHY THIS EXISTS.  The `sinks` block of tools/model_report.c cannot be trusted
 * and cannot be repaired in place.  Three defects:
 *
 * 1. STALE ALGEBRAIC STATE.  It rebuilds each step as `Input s = in;` where `in`
 *    is the Input that simulate_days() has already mutated to the FINAL step of
 *    the run.  Only sucrose, starch and photosynthesis are refreshed from the
 *    recorded trajectory.  Every other algebraic quantity --
 *    starch_degradation_rate, uptake_cost, transport_cost, nitrogen_uptake,
 *    sucrose_loading_frequency -- stays frozen at its end-of-run value.
 *    starch_degradation_rate carries a factor (1 - light), so freezing it
 *    applies a night-only flux across the light period too.
 *
 * 2. STATE THE TRAJECTORY DOES NOT CARRY.  SimulationResult holds only sucrose,
 *    starch, photosynthesis and total biomass.  N and P are integrated but not
 *    recorded, so no post-hoc reconstruction can see them.  This probe therefore
 *    does NOT reconstruct: it recomputes every flux inside the observer, from
 *    the live solver state vector, which does carry X_N and X_P.  If uptake
 *    still reads zero here, that is a property of the model, not of the tool --
 *    and the uptake diagnostics below are written out so the two can be told
 *    apart.
 *
 * 3. A DOUBLE COUNT WITH A CONSERVATION PROOF.  Over one diel cycle in periodic
 *    steady state, integrate the two carbohydrate ODEs of dynamics.c:
 *
 *      starch :  0 = int k_s A dt  -  int D dt
 *      sucrose:  0 = int (1-k_s) A dt + int D dt
 *                    - int [uc + tc + r_m S + lambda_L nes] dt
 *
 *    Substituting the first into the second gives, exactly,
 *
 *      int A dt = int [uc + tc + r_m S + lambda_L nes] dt
 *
 *    THE FOUR SUCROSE SINKS ALONE ACCOUNT FOR ALL GROSS ASSIMILATION.  Starch is
 *    not a fifth sink at the daily scale; it is an intra-day buffer whose net
 *    flux over a closed cycle is zero, and its synthesis is already inside the
 *    four sucrose sinks by way of the night degradation term.  Reporting starch
 *    as a fifth share and summing is therefore guaranteed to overshoot 1.0,
 *    which is the ~1.9x that has been reported.
 *
 * WHAT THIS TOOL REPORTS INSTEAD.
 *    - the four sucrose sinks as shares of GROSS assimilation; these must sum to
 *      1.0, and the residual is written out as closure_residual so the
 *      accounting is falsifiable rather than assumed
 *    - starch GROSS synthesis and GROSS degradation as separate shares.  Gross
 *      synthesis is what Sulpice 2014 mmc3 row 24 "% assimilation rate to starch
 *      synth" measures, so it is the only starch number comparable to the paper
 *    - starch_net_share, which must be ~0 if the cycle is truly periodic; a
 *      large value means the run has not converged and every share is void
 *    - uptake diagnostics (N pool, N uptake, uptake_cost) so a zero uptake share
 *      can be attributed
 */

#include "config/config.h"
#include "model/input.h"
#include "model/model.h"
#include "model/state.h"
#include "simulation/dynamics.h"
#include "simulation/observers.h"
#include "simulation/simulation.h"

#include <solver.h>
#ifndef max_align_t
#define max_align_t long double
#define BP_MAX_ALIGN_T_SHIM
#endif
#include <framework/default_runtime.h>
#include <framework/model_buffers.h>
#include <framework/model_runner.h>
#ifdef BP_MAX_ALIGN_T_SHIM
#undef max_align_t
#undef BP_MAX_ALIGN_T_SHIM
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>


/* step budget, same mechanism as the ported one in simulation.c -- the probe
 * builds its own solver runtime, so it needs its own registration. */
#define BP_MAX_SOLVER_STEPS 10000000u
static unsigned long bp_step_budget(void)
{
    const char* s = getenv("PLANTMODEL_MAX_SOLVER_STEPS");
    if (s && *s) { char* e=NULL; unsigned long v=strtoul(s,&e,10); if (e!=s && v>0ul) return v; }
    return (unsigned long)BP_MAX_SOLVER_STEPS;
}
typedef struct { unsigned long steps, budget; int exhausted; } BpBudget;
static SolverPluginResult bp_on_attempt(void* vctx, const SolverEvent* ev)
{
    (void)ev; BpBudget* b=(BpBudget*)vctx;
    if (!b) return SOLVER_PLUGIN_CONTINUE;
    if (++b->steps >= b->budget) { b->exhausted = 1; return SOLVER_PLUGIN_STOP; }
    return SOLVER_PLUGIN_CONTINUE;
}

#define NPP 4
static const int    PP[NPP]        = { 4, 6, 8, 12 };
static const double LAMBDA_C[NPP]  = { 0.82, 0.79, 0.67, 0.62 };
static const double LAMBDA_SNI[NPP]= { 0.16, 0.15, 0.13, 0.08 };
static const double LAMBDA_SB[NPP] = { 0.003843, 0.004433, 0.005246, 0.005850 };

/* Measured budget, for comparison only -- nothing here is fitted to it.
 * Starch gross synthesis share: Sulpice 2014 mmc3 "Growth components" row 24. */
static const double STARCH_SYNTH_MEASURED[NPP] = { 0.4861, 0.5129, 0.4143, 0.3272 };

/* ------------------------------------------------------- per-step recording */

#define MAXS 40000
typedef struct {
    real_t next_t, dt;
    int    n;
    double dt_used;
    PlantCtx* ctx;
    void (*unpack)(Input*, const real_t*);
    void (*algeb)(Input*, real_t, real_t);
    double t[MAXS], light[MAXS], A[MAXS], ks[MAXS], D[MAXS];
    double uc[MAXS], tc[MAXS], rm[MAXS], ld[MAXS];
    double suc[MAXS], sta[MAXS];
    double npool[MAXS], nupt[MAXS], ppool[MAXS], pupt[MAXS];
} BudgetObs;

static void budget_observer(real_t t, const real_t* x, size_t n, void* vctx)
{
    (void)n;
    BudgetObs* o = (BudgetObs*)vctx;
    if (t + REAL(1e-15) < o->next_t) return;
    while (o->n < MAXS && o->next_t <= t + REAL(1e-15)) {
        Input s = *(o->ctx->base);
        o->unpack(&s, x);
        /* recompute EVERY algebraic quantity at this hour from the live state --
         * the step model_report's reconstruction skips */
        o->algeb(&s, fmod(o->next_t, REAL(24.0)), o->ctx->starch_night_start);

        real_t nes = night_efficiency_starch(s.carbohydrates.sucrose,
                                            s.carbohydrates.max_sucrose,
                                            s.carbohydrates.lambda_g,
                                            s.photo.light,
                                            s.carbohydrates.starch_partition_coeff,
                                            s.carbohydrates.sink_activation_half,
                                            s.carbohydrates.sink_activation_hill);
        int i = o->n;
        o->t[i]     = (double)o->next_t;
        o->light[i] = (double)s.photo.light;
        o->A[i]     = (double)s.photo.photosynthesis;
        o->ks[i]    = (double)s.carbohydrates.starch_partition_coeff;
        o->D[i]     = (double)s.carbohydrates.starch_degradation_rate;
        o->uc[i]    = (double)s.carbohydrates.uptake_cost;
        o->tc[i]    = (double)s.carbohydrates.transport_cost;
        o->rm[i]    = (double)s.carbohydrates.respiration_frequency
                      * (double)s.carbohydrates.sucrose;
        o->ld[i]    = (double)s.carbohydrates.sucrose_loading_frequency * (double)nes;
        o->suc[i]   = (double)s.carbohydrates.sucrose;
        o->sta[i]   = (double)s.carbohydrates.starch;
        o->npool[i] = (double)s.nutrients.nitrogen;
        o->nupt[i]  = (double)s.nutrients.nitrogen_uptake;
        o->ppool[i] = (double)s.nutrients.phosphorus;
        o->pupt[i]  = (double)s.nutrients.phosphorus_uptake;
        o->n++;
        o->next_t += o->dt;
    }
}

/* ------------------------------------------------------------------- runner */

typedef struct { double kappa; int chamber, days; } Variant;

static int setup(int i, const Variant* v, double lsb, SimulationConfig* cfg, Input* in)
{
    *cfg = simulation_default_config();
    *in  = generate_input();
    if (config_load_file("config/default.config", cfg, in) != 0) return -1;
    cfg->days = v->days; cfg->write_csv = 0; cfg->gui_enabled = 0;
    in->photo.photo_mode                 = (int)PHOTO_MODE_FARQUHAR_JARVIS;
    in->photo.photoperiod                = (real_t)PP[i];
    in->photo.feedback_on_photosynthesis = (real_t)LAMBDA_C[i];
    in->carbohydrates.lambda_sni         = (real_t)LAMBDA_SNI[i];
    in->growth.lambda_sb                 = (real_t)lsb;
    /* kappa enters through SLA, exactly as model_report does it, so the two
     * tools are measuring the same variant. */
    in->gas_exchange.specific_leaf_area  = (real_t)(0.025 * (v->kappa / 1.10));
    if (v->chamber) {
        in->photo.max_PAR                 = REAL(160.0);
        in->gas_exchange.leaf_temperature = REAL(21.0);
    }
    return 0;
}

static int run_with_budget(int i, const Variant* v, double lsb, BudgetObs* obs)
{
    SimulationConfig cfg; Input in;
    if (setup(i, v, lsb, &cfg, &in) != 0) return -1;

    PlantCtx ctx = { .base = &in, .starch_night_start = in.carbohydrates.starch, .was_light = 0 };
    NightStartObserverCtx nctx = { .plant_ctx = &ctx, .starch_index = X_STARCH };

    memset(obs, 0, sizeof(*obs));
    obs->next_t = REAL(0.0); obs->dt = cfg.dt_hours; obs->ctx = &ctx;
    obs->unpack = simulation_unpack_state; obs->algeb = simulation_compute_algebraic;
    obs->dt_used = (double)cfg.dt_hours;

    VectorObserver obs_list[3] = { nan_guard_observer, night_start_observer, budget_observer };
    void* ctx_list[3] = { NULL, &nctx, obs };
    ObserverChain chain = { .obs = obs_list, .ctx = ctx_list, .count = 3 };

    SolverStats stats; solver_stats_init(&stats);
    SolverPluginManager plugins; solver_plugin_manager_init(&plugins);
    SolverDefaultRuntimeConfig rc = { .observer = observer_chain_cb,
                                      .observer_ctx = &chain,
                                      .stats = &stats, .logger = NULL };
    if (solver_default_runtime_init(&plugins, &rc) != 0) plugins.count = 0;

    static BpBudget bud;
    bud = (BpBudget){ .steps = 0, .budget = bp_step_budget(), .exhausted = 0 };
    SolverPlugin bp = { .name = "budget_probe_step_budget", .plugin_ctx = &bud,
                        .on_run_start = NULL, .on_step_attempt = bp_on_attempt,
                        .on_step_accepted = NULL, .on_step_rejected = NULL,
                        .on_run_finish = NULL };
    solver_plugin_manager_add(&plugins, &bp);

    const ModelInterface* model = simulation_plant_model_interface();
    ModelInterface rm_model = *model; rm_model.model_ctx = &ctx;
    ModelBuffers buf; memset(&buf, 0, sizeof(buf));
    SolverRunConfig run = { .type = cfg.solver_type, .t0 = REAL(0.0),
                            .t_end = (real_t)cfg.days * REAL(24.0),
                            .h_init = cfg.dt_hours, .tol = cfg.solver_tolerance,
                            .plugins = &plugins };
    if (model_buffers_init(&buf, rm_model.descriptor) != MODEL_STATUS_OK) return -1;
    simulation_pack_state(buf.state, &in);
    ModelStatus st = model_run(&rm_model, &buf, &run);
    model_buffers_free(&buf);
    if (st != MODEL_STATUS_OK) return -1;
    if (bud.exhausted) return -1;   /* stalled adaptive step = divergent */
    for (int k = 0; k < obs->n; ++k)
        if (!isfinite(obs->A[k]) || !isfinite(obs->suc[k]) || !isfinite(obs->sta[k])) return -1;
    return 0;
}

static double trap(const double* y, int a, int b, double dt)
{
    double s = 0.0;
    for (int k = a; k < b - 1; ++k) s += 0.5 * (y[k] + y[k + 1]) * dt;
    return s;
}

static void ensure_dir(const char* p)
{
    char b[512]; snprintf(b, sizeof(b), "%s", p);
    for (char* q = b + 1; *q; ++q) if (*q == '/') { *q = 0; mkdir(b, 0755); *q = '/'; }
}

int main(int argc, char** argv)
{
    Variant v = { .kappa = 1.10, .chamber = 0, .days = 30 };
    const char* out = "results/attribution/budget.csv";
    const char* label = "budget";
    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "--out") && i + 1 < argc)        out = argv[++i];
        else if (!strcmp(argv[i], "--label") && i + 1 < argc) label = argv[++i];
        else if (!strcmp(argv[i], "--kappa") && i + 1 < argc) v.kappa = atof(argv[++i]);
        else if (!strcmp(argv[i], "--days") && i + 1 < argc)  v.days = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--chamber"))               v.chamber = 1;
    }
    ensure_dir(out);
    FILE* fp = fopen(out, "w");
    if (!fp) { fprintf(stderr, "cannot open %s\n", out); return 1; }
    fprintf(fp, "label,kappa,chamber,days,photoperiod,metric,value,measured_lo,measured_hi\n");

    static BudgetObs obs;
    for (int i = 0; i < NPP; ++i) {
        if (run_with_budget(i, &v, LAMBDA_SB[i], &obs) != 0) {
            fprintf(fp, "%s,%.4f,%d,%d,%d,diverged,1,,\n",
                    label, v.kappa, v.chamber, v.days, PP[i]);
            continue;
        }
        /* dt comes from the config (0.1 h), NOT assumed hourly: an assumed dt of
         * 1.0 indexes a window ten times too short, landing wholly in the night
         * and reading as zero assimilation. */
        double dt = obs.dt_used;
        int spd = (int)lround(24.0 / dt);
        int c0 = (v.days - 1) * spd, c1 = v.days * spd + 1;
        if (c1 > obs.n) c1 = obs.n;
        int dusk = c0 + (int)lround(PP[i] / dt);
        if (dusk >= obs.n) dusk = obs.n - 1;
        int d1 = dusk + 1; if (d1 > obs.n) d1 = obs.n;

        double A_cyc = trap(obs.A, c0, c1, dt);
        double A_lit = trap(obs.A, c0, d1, dt);

        /* the four sucrose sinks -- these are the whole budget (header proof) */
        double uc = trap(obs.uc, c0, c1, dt);
        double tc = trap(obs.tc, c0, c1, dt);
        double rm = trap(obs.rm, c0, c1, dt);
        double ld = trap(obs.ld, c0, c1, dt);

        /* starch: gross synthesis and gross degradation, kept separate */
        static double gs[MAXS];
        for (int k = c0; k < c1; ++k) gs[k] = obs.ks[k] * obs.A[k];
        double S_syn = trap(gs, c0, c1, dt);
        double S_deg = trap(obs.D, c0, c1, dt);

        double sum4 = uc + tc + rm + ld;
        #define W(m, val, lo, hi) fprintf(fp, "%s,%.4f,%d,%d,%d,%s,%.6f,%s,%s\n", \
            label, v.kappa, v.chamber, v.days, PP[i], m, (val), lo, hi)
        W("assim_gross_cycle", A_cyc, "", "");
        W("assim_gross_light", A_lit, "", "");
        W("assim_mean_light_rate", A_lit / (double)PP[i], "", "");
        W("share_uptake",     uc / A_cyc, "", "");
        W("share_transport",  tc / A_cyc, "", "");
        W("share_respiration",rm / A_cyc, "0.10", "0.18");
        W("share_phloem_load",ld / A_cyc, "", "");
        W("sum_four_sinks",   sum4 / A_cyc, "1.0", "1.0");
        W("closure_residual", (sum4 - A_cyc) / A_cyc, "0", "0");
        fprintf(fp, "%s,%.4f,%d,%d,%d,starch_gross_synth_share,%.6f,%.4f,%.4f\n",
                label, v.kappa, v.chamber, v.days, PP[i], S_syn / A_cyc,
                STARCH_SYNTH_MEASURED[i], STARCH_SYNTH_MEASURED[i]);
        W("starch_gross_deg_share", S_deg / A_cyc, "", "");
        W("starch_net_share",       (S_syn - S_deg) / A_cyc, "0", "0");
        W("dusk_starch",  obs.sta[dusk], "", "");
        W("dusk_sucrose", obs.suc[dusk], "", "");
        /* uptake diagnostics: distinguish "model does not take up N" from
         * "the tool cannot see N" */
        W("n_pool_dusk",      obs.npool[dusk], "", "");
        W("n_uptake_dusk",    obs.nupt[dusk], "", "");
        W("n_uptake_cycle",   trap(obs.nupt, c0, c1, dt), "", "");
        W("p_pool_dusk",      obs.ppool[dusk], "", "");
        W("p_uptake_cycle",   trap(obs.pupt, c0, c1, dt), "", "");
        W("uptake_cost_dusk", obs.uc[dusk], "", "");
        #undef W
    }
    fclose(fp);
    printf("wrote %s (label=%s kappa=%.4f chamber=%d)\n", out, label, v.kappa, v.chamber);
    return 0;
}

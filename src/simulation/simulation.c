#include "export/export.h"
#include "gui/plot.h"
#include "model/input.h"
#include "model/model.h"
#include "simulation/dynamics.h"
#include "simulation/observers.h"
#include "simulation/simulation.h"
#include <solver.h>
#include <logger.h>
#include <observers.h>
#include <pthread.h>
#include <stddef.h>
#ifndef max_align_t
#define max_align_t long double
#define PLANTMODEL_MAX_ALIGN_T_SHIM
#endif
#include <framework/default_runtime.h>
#include <framework/model_buffers.h>
#include <framework/model_runner.h>
#ifdef PLANTMODEL_MAX_ALIGN_T_SHIM
#undef max_align_t
#undef PLANTMODEL_MAX_ALIGN_T_SHIM
#endif

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

static const real_t DEFAULT_SIMULATION_DT = REAL(1);

/*
 * Step budget for the adaptive integrator.
 *
 * Failure mode this guards against: for some parameter combinations the
 * adaptive step collapses without the state ever becoming non-finite, so the
 * integrator neither advances nor reports an error and simulate_days() never
 * returns.  Observed at the 12 h photoperiod for lambda_sdi above roughly
 * 0.516 (threshold bisected to [0.5162, 0.5165] at lambda_sb = 6.66e-3): the
 * run issues tens of millions of gas-exchange solves per wall-clock second,
 * every one of them convergent, while t stands still.  The NaN guard cannot
 * see this because nothing is NaN, and neither can a divergence test on the
 * outputs, because no outputs are ever produced.
 *
 * A step budget converts the hang into a classified divergence: the plugin
 * returns SOLVER_PLUGIN_STOP, the run unwinds normally, and the caller sees
 * result->diverged, exactly as it would for a non-finite trajectory.  Batch
 * callers (sensitivity, cross-validation, sweeps) can then record the
 * parameter combination as inadmissible and carry on, which is what they
 * already do for divergence.
 *
 * The budget is generous relative to legitimate work: production runs use
 * dt = 0.1 h over at most 35 days, i.e. 8400 output samples, and the adaptive
 * integrator takes a low multiple of that in internal steps.  Ten million
 * steps is more than three orders of magnitude of headroom, so no admissible
 * run can reach it, while a stalled run reaches it in seconds.
 */
#define SIMULATION_MAX_SOLVER_STEPS 10000000u

/*
 * Optional override, read from the environment as PLANTMODEL_MAX_SOLVER_STEPS.
 *
 * The default budget is sized so that no admissible run can reach it, which
 * means a stalled run burns ten million steps before it is classified.  That
 * is the right trade for a production run, where a false divergence would be
 * far worse than a slow one, but it is the wrong trade inside a parameter
 * search: a Nelder-Mead simplex that wanders into the stiff region spends
 * minutes per stalled evaluation, and the fit becomes dominated by points it
 * is going to reject anyway.
 *
 * The override lets a fitting harness lower the budget to something a few
 * times the legitimate step count, so stalled candidates are rejected in
 * milliseconds.  It cannot change the outcome for an admissible run, because
 * an admissible run never approaches the budget -- it only changes how long
 * the model takes to say no.  Default behaviour is unchanged when the variable
 * is unset, so production and validation runs are unaffected.
 */
static unsigned long simulation_step_budget(void)
{
    const char* s = getenv("PLANTMODEL_MAX_SOLVER_STEPS");
    if (s && *s) {
        char* end = NULL;
        unsigned long v = strtoul(s, &end, 10);
        if (end != s && v > 0ul) return v;
    }
    return (unsigned long)SIMULATION_MAX_SOLVER_STEPS;
}

typedef struct {
    unsigned long steps;
    unsigned long budget;
    int           exhausted;
} StepBudgetCtx;

static SolverPluginResult step_budget_on_attempt(void* vctx, const SolverEvent* event)
{
    (void)event;
    StepBudgetCtx* b = (StepBudgetCtx*)vctx;
    if (!b)
        return SOLVER_PLUGIN_CONTINUE;
    if (++b->steps >= b->budget) {
        if (!b->exhausted) {
            b->exhausted = 1;
            fprintf(stderr,
                "warning: solver step budget (%lu) exhausted at t=%.6f h; the "
                "adaptive step has collapsed without the state becoming "
                "non-finite. Stopping and classifying this run as divergent.\n",
                b->budget, (double)event->t_start);
        }
        return SOLVER_PLUGIN_STOP;
    }
    return SOLVER_PLUGIN_CONTINUE;
}

typedef struct {
    char filename[64];
    CSVObserverCtx csv_ctx;
    SimObsCtx sample_ctx;
    NightStartObserverCtx night_ctx;
    NanGuardCtx nan_ctx;
    DomainGuardCtx domain_ctx;
    StepBudgetCtx budget_ctx;
    VectorObserver obs_list[5];
    void* ctx_list[5];
    ObserverChain chain;
    SolverStats stats;
    SolverPluginManager plugins;
    SolverLogger logger;
    SolverLogger* logger_ptr;
    GuiOutput gui_output;
    GuiThreadArgs gui_args;
    pthread_t gui_thread;
    int gui_started;
} SimulationRuntime;

static void simulation_sync_buffers_from_input(ModelBuffers* buffers, const Input* input)
{
    simulation_pack_state(buffers->state, input);
}

static void simulation_sync_input_from_buffers(Input* input, const ModelBuffers* buffers)
{
    simulation_unpack_state(input, buffers->state);
}

static SimulationConfig simulation_normalize_config(const SimulationConfig* config)
{
    SimulationConfig cfg = config ? *config : simulation_default_config();
    if (cfg.days <= 0)
        cfg.days = 1;
    if (cfg.dt_hours <= REAL(0.0))
        cfg.dt_hours = DEFAULT_SIMULATION_DT;
    return cfg;
}

static void simulation_result_reset(SimulationResult* result)
{
    result->sucrose         = NULL;
    result->starch          = NULL;
    result->ph              = NULL;
    result->total_biomass       = NULL;
    result->sucrose_vacuole = NULL;
    result->sucrose_total   = NULL;
    result->filled_steps    = 0;
    result->total_steps     = 0;
    result->photoperiod     = REAL(0.0);
    result->days            = 0;
    result->diverged        = 0;
    result->diverged_t      = 0.0;
    result->diverged_index  = -1;
    result->stalled           = 0;
    result->left_domain       = 0;
    result->domain_exit_t     = 0.0;
    result->domain_exit_index = -1;
    result->domain_min        = 0.0;
}

static void simulation_result_init(SimulationResult* result, int days, real_t photoperiod, real_t dt_hours)
{
    simulation_result_reset(result);
    result->days = days;
    result->photoperiod = photoperiod;
    result->total_steps = (int)(days * REAL(24.0) / dt_hours);
    result->sucrose         = calloc((size_t)result->total_steps, sizeof(real_t));
    result->starch          = calloc((size_t)result->total_steps, sizeof(real_t));
    result->ph              = calloc((size_t)result->total_steps, sizeof(real_t));
    result->total_biomass       = calloc((size_t)result->total_steps, sizeof(real_t));
    result->sucrose_vacuole = calloc((size_t)result->total_steps, sizeof(real_t));
    result->sucrose_total   = calloc((size_t)result->total_steps, sizeof(real_t));
}

static void simulation_runtime_init_csv(SimulationRuntime* runtime, const SimulationConfig* cfg, real_t photoperiod)
{
    snprintf(runtime->filename, sizeof(runtime->filename), "results/simulation/plant_%dh.csv", (int)photoperiod);
    runtime->csv_ctx.f = cfg->write_csv ? fopen(runtime->filename, "w") : NULL;
    /* TODO: Surface CSV open failures to caller/tests (currently silently continues without output file). */

    if (runtime->csv_ctx.f) {
        fprintf(runtime->csv_ctx.f,
            "t,light,photosynthesis,"
            "leaf_biomass,root_biomass,total_biomass,"
            "sucrose,starch,"
            "nitrogen,phosphorus,"
            "nitrogen_uptake,phosphorus_uptake,"
            "uptake_cost,transport_cost,"
            "starch_degradation_rate,night_efficiency\n");

        fflush(runtime->csv_ctx.f);
    }
}

static void simulation_runtime_init_logging(SimulationRuntime* runtime, const SimulationConfig* cfg)
{
    solver_stats_init(&runtime->stats);
    solver_plugin_manager_init(&runtime->plugins);
    runtime->logger_ptr = NULL;

    if (cfg->write_csv) {
        if (solver_logger_init(&runtime->logger, "steps.csv", "summary.csv") == 0) {
            runtime->logger_ptr = &runtime->logger;
        }
    }

    SolverDefaultRuntimeConfig runtime_cfg = {
        .observer = observer_chain_cb,
        .observer_ctx = &runtime->chain,
        .stats = &runtime->stats,
        .logger = runtime->logger_ptr
    };
    if (solver_default_runtime_init(&runtime->plugins, &runtime_cfg) != 0) {
        /* Fallback: run without runtime plugins if default runtime setup fails. */
        runtime->plugins.count = 0;
    }

    /* Step budget: turns a collapsed adaptive step into a classified
     * divergence instead of a hang.  Registered last so it runs after the
     * default runtime's own hooks, and pointing at runtime-owned context so
     * the counter is per-run rather than per-process. */
    runtime->budget_ctx = (StepBudgetCtx){
        .steps = 0, .budget = simulation_step_budget(), .exhausted = 0
    };
    SolverPlugin budget_plugin = {
        .name             = "plantmodel_step_budget",
        .plugin_ctx       = &runtime->budget_ctx,
        .on_run_start     = NULL,
        .on_step_attempt  = step_budget_on_attempt,
        .on_step_accepted = NULL,
        .on_step_rejected = NULL,
        .on_run_finish    = NULL
    };
    solver_plugin_manager_add(&runtime->plugins, &budget_plugin);
}

static void simulation_runtime_init_gui(
    SimulationRuntime* runtime,
    const SimulationConfig* cfg,
    SimulationResult* result)
{
    runtime->gui_output = (GuiOutput){
        .sucrose = result->sucrose,
        .starch = result->starch,
        .photosynthesis = result->ph,
        .total_biomass = result->total_biomass,
        .filled_steps = &result->filled_steps,
        .total_steps = result->total_steps,
        .dt_hours = cfg->dt_hours,
        .photoperiod = result->photoperiod,
        .days = result->days
    };
    runtime->gui_args = (GuiThreadArgs){
        .output = &runtime->gui_output,
        .title = "Plant Simulation",
        .target_fps = 60,
        .summary_csv_path = "summary.csv"
    };

    runtime->gui_started = 0;
    if (cfg->gui_enabled) {
        runtime->gui_started = (pthread_create(&runtime->gui_thread, NULL, gui_main_thread, &runtime->gui_args) == 0);
    }
    /* TODO: Surface GUI thread startup failures so caller can react (headless fallback/logging). */
}

static void simulation_runtime_cleanup(SimulationRuntime* runtime)
{
    if (runtime->csv_ctx.f)
        fclose(runtime->csv_ctx.f);

    if (runtime->logger_ptr) {
        solver_logger_close(&runtime->logger);
    }

    if (runtime->gui_started)
        pthread_join(runtime->gui_thread, NULL);
}

static void simulation_prepare_runtime(
    SimulationRuntime* runtime,
    const SimulationConfig* cfg,
    Input* input,
    SimulationResult* result,
    PlantCtx* rhs_ctx)
{
    *runtime = (SimulationRuntime){0};

    runtime->csv_ctx = (CSVObserverCtx){
        .f = NULL,
        .plant_ctx = rhs_ctx,
        .unpack_state = simulation_unpack_state,
        .compute_algebraic = simulation_compute_algebraic
    };
    runtime->sample_ctx = (SimObsCtx){
        .next_t = REAL(0.0),
        .dt_sample = cfg->dt_hours,
        .step = 0,
        .max_steps = result->total_steps,
        .rhs_ctx = rhs_ctx,
        .result = result,
        .unpack_state = simulation_unpack_state,
        .compute_algebraic = simulation_compute_algebraic
    };
    runtime->night_ctx = (NightStartObserverCtx){
        .plant_ctx = rhs_ctx,
        .starch_index = X_STARCH
    };
    runtime->nan_ctx = (NanGuardCtx){
        .tripped = 0, .t_first = 0.0, .index_first = 0, .n_events = 0, .warned = 0
    };
    runtime->domain_ctx = (DomainGuardCtx){
        .tripped = 0, .t_first = 0.0, .index_first = 0,
        .min_value = 0.0, .index_min = 0, .n_events = 0
    };
    runtime->obs_list[0] = nan_guard_observer;
    runtime->obs_list[1] = domain_guard_observer;
    runtime->obs_list[2] = night_start_observer;
    runtime->obs_list[3] = sim_sampling_observer;
    runtime->obs_list[4] = cfg->write_csv ? csv_observer : NULL;
    runtime->ctx_list[0] = &runtime->nan_ctx;
    runtime->ctx_list[1] = &runtime->domain_ctx;
    runtime->ctx_list[2] = &runtime->night_ctx;
    runtime->ctx_list[3] = &runtime->sample_ctx;
    runtime->ctx_list[4] = cfg->write_csv ? (void*)&runtime->csv_ctx : NULL;
    runtime->chain = (ObserverChain){
        .obs = runtime->obs_list,
        .ctx = runtime->ctx_list,
        .count = 5
    };

    simulation_runtime_init_csv(runtime, cfg, input->photo.photoperiod);
    simulation_runtime_init_logging(runtime, cfg);
    simulation_runtime_init_gui(runtime, cfg, result);
}

static ModelStatus simulation_run_model(
    const SimulationConfig* cfg,
    Input* input,
    ModelInterface* runtime_model,
    ModelBuffers* buffers,
    SimulationRuntime* runtime)
{
    SolverRunConfig run_cfg = {
        .type = cfg->solver_type,
        .t0 = REAL(0.0),
        .t_end = (real_t)cfg->days * REAL(24.0),
        .h_init = cfg->dt_hours,
        .tol = cfg->solver_tolerance,
        .plugins = &runtime->plugins
    };
    ModelStatus status = model_buffers_init(buffers, runtime_model->descriptor);
    if (status != MODEL_STATUS_OK) {
        fprintf(stderr, "model_buffers_init failed: %s\n", model_status_string(status));
        return status;
    }

    simulation_sync_buffers_from_input(buffers, input);
    status = model_run(runtime_model, buffers, &run_cfg);
    if (status != MODEL_STATUS_OK) {
        fprintf(stderr, "model_run failed: %s\n", model_status_string(status));
        return status;
    }

    simulation_sync_input_from_buffers(input, buffers);
    return MODEL_STATUS_OK;
}

void simulate_step(real_t current_time, int step, Input* input, int days)
{
    (void)step;
    /* TODO: Remove or wire this function into production flow; it is currently unused by the main simulation path. */
    real_t hour_of_day = fmod(current_time, REAL(24.0));
    update_light_conditions(input, hour_of_day);

    static bool trigger = false;
    static real_t starch_night_start = REAL(0.0);

    if (input->photo.light == REAL(0.0)) {
        if (!trigger) {
            starch_night_start = input->carbohydrates.starch;
        }
        trigger = true;
    } else {
        trigger = false;
    }

    const ModelInterface* model = simulation_plant_model_interface();
    ModelBuffers buffers = {0};
    ModelStatus status;

    PlantCtx ctx = {
        .base = input,
        .starch_night_start = starch_night_start
    };
    ModelInterface runtime_model = *model;
    runtime_model.model_ctx = &ctx;

    SolverRunConfig run_cfg = {
        .type = SOLVER_ODE45,
        .t0 = current_time,
        .t_end = current_time + DEFAULT_SIMULATION_DT,
        .h_init = DEFAULT_SIMULATION_DT,
        .tol = REAL(1e-8),
        .plugins = NULL
    };

    status = model_buffers_init(&buffers, runtime_model.descriptor);
    if (status != MODEL_STATUS_OK) {
        fprintf(stderr, "model_buffers_init failed: %s\n", model_status_string(status));
        return;
    }

    simulation_sync_buffers_from_input(&buffers, input);
    status = model_run(&runtime_model, &buffers, &run_cfg);
    if (status != MODEL_STATUS_OK) {
        fprintf(stderr, "model_run failed: %s\n", model_status_string(status));
        model_buffers_free(&buffers);
        return;
    }

    simulation_sync_input_from_buffers(input, &buffers);
    model_buffers_free(&buffers);

    simulation_compute_algebraic(
        input,
        fmod(current_time + DEFAULT_SIMULATION_DT, REAL(24.0)),
        starch_night_start);
    log_simulation_step(days, current_time, input);
}

SimulationConfig simulation_default_config(void)
{
    SimulationConfig cfg = {
        .days = 30,
        .dt_hours = DEFAULT_SIMULATION_DT,
        .solver_tolerance = REAL(1e-8),
        .solver_type = SOLVER_ODE45,
        .gui_enabled = 1,
        .write_csv = 1
    };
    return cfg;
}

void simulate_days(const SimulationConfig* config, Input* input, SimulationResult* result)
{
    SimulationConfig cfg = simulation_normalize_config(config);
    simulation_result_init(result, cfg.days, input->photo.photoperiod, cfg.dt_hours);
    /* TODO: Validate allocation results and return an explicit error code instead of assuming success. */

    input->growth.lambda_sb = lambda_sb_f(
        input->growth.lambda_sb,
        input->nutrients.nitrogen_soil_content,
        input->nutrients.phosphorus_soil_content);

    const ModelInterface* model = simulation_plant_model_interface();
    ModelInterface runtime_model = *model;
    ModelBuffers buffers = {0};
    ModelStatus status = MODEL_STATUS_OK;

    PlantCtx rhs_ctx = {
        .base = input,
        .starch_night_start = input->carbohydrates.starch
    };
    runtime_model.model_ctx = &rhs_ctx;

    Input tmp = *input;
    update_light_conditions(&tmp, REAL(0.0));
    rhs_ctx.was_light = (tmp.photo.light != REAL(0.0));

    SimulationRuntime runtime;
    simulation_prepare_runtime(&runtime, &cfg, input, result, &rhs_ctx);
    status = simulation_run_model(&cfg, input, &runtime_model, &buffers, &runtime);
    if (status == MODEL_STATUS_OK)
        model_buffers_free(&buffers);

    if (status != MODEL_STATUS_OK && buffers.state)
        model_buffers_free(&buffers);

    /* Surface the divergence record so callers can classify the run.
     * An exhausted step budget counts as divergence: the trajectory never
     * reached t_end, so the result is not a solution of the system even though
     * every sampled state stayed finite. */
    result->diverged       = runtime.nan_ctx.tripped || runtime.budget_ctx.exhausted;
    result->diverged_t     = runtime.nan_ctx.t_first;
    result->diverged_index = runtime.nan_ctx.tripped
                                 ? (int)runtime.nan_ctx.index_first
                                 : -1;

    /* Break `diverged` down into the reason, without changing what `diverged`
     * itself means to existing callers.  `stalled` is the budget case with a
     * finite state throughout: the partial trajectory is a valid solution and
     * a caller that only needs behaviour before diverged_t may use it.
     * `left_domain` is orthogonal and is reported even for runs that complete,
     * because a completed run whose pools went negative is finite but outside
     * the model's domain of validity. */
    result->stalled           = runtime.budget_ctx.exhausted && !runtime.nan_ctx.tripped;
    result->left_domain       = runtime.domain_ctx.tripped;
    result->domain_exit_t     = runtime.domain_ctx.t_first;
    result->domain_exit_index = runtime.domain_ctx.tripped
                                    ? (int)runtime.domain_ctx.index_first
                                    : -1;
    result->domain_min        = runtime.domain_ctx.min_value;

    simulation_runtime_cleanup(&runtime);
}

void simulation_result_free(SimulationResult* result)
{
    if (!result)
        return;
    free(result->sucrose);
    free(result->starch);
    free(result->ph);
    free(result->total_biomass);
    free(result->sucrose_vacuole);
    free(result->sucrose_total);
    result->sucrose         = NULL;
    result->starch          = NULL;
    result->ph              = NULL;
    result->total_biomass       = NULL;
    result->sucrose_vacuole = NULL;
    result->sucrose_total   = NULL;
}

real_t compute_rgr(real_t FW_start, real_t FW_end, real_t duration_hours)
{
    return (RLOG(FW_end) - RLOG(FW_start)) / duration_hours;
}

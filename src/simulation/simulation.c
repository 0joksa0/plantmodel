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

typedef struct {
    char filename[64];
    CSVObserverCtx csv_ctx;
    SimObsCtx sample_ctx;
    NightStartObserverCtx night_ctx;
    VectorObserver obs_list[4];
    void* ctx_list[4];
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
    result->sucrose = NULL;
    result->starch = NULL;
    result->ph = NULL;
    result->partition = NULL;
    result->filled_steps = 0;
    result->total_steps = 0;
    result->photoperiod = REAL(0.0);
    result->days = 0;
}

static void simulation_result_init(SimulationResult* result, int days, real_t photoperiod, real_t dt_hours)
{
    simulation_result_reset(result);
    result->days = days;
    result->photoperiod = photoperiod;
    result->total_steps = (int)(days * REAL(24.0) / dt_hours);
    result->sucrose = calloc((size_t)result->total_steps, sizeof(real_t));
    result->starch = calloc((size_t)result->total_steps, sizeof(real_t));
    result->ph = calloc((size_t)result->total_steps, sizeof(real_t));
    result->partition = calloc((size_t)result->total_steps, sizeof(real_t));
}

static void simulation_runtime_init_csv(SimulationRuntime* runtime, const SimulationConfig* cfg, real_t photoperiod)
{
    snprintf(runtime->filename, sizeof(runtime->filename), "plant_%dh.csv", (int)photoperiod);
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
        .total_biomass = result->partition,
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
    runtime->obs_list[0] = nan_guard_observer;
    runtime->obs_list[1] = night_start_observer;
    runtime->obs_list[2] = sim_sampling_observer;
    runtime->obs_list[3] = cfg->write_csv ? csv_observer : NULL;
    runtime->ctx_list[0] = NULL;
    runtime->ctx_list[1] = &runtime->night_ctx;
    runtime->ctx_list[2] = &runtime->sample_ctx;
    runtime->ctx_list[3] = cfg->write_csv ? (void*)&runtime->csv_ctx : NULL;
    runtime->chain = (ObserverChain){
        .obs = runtime->obs_list,
        .ctx = runtime->ctx_list,
        .count = 4
    };

    simulation_runtime_init_csv(runtime, cfg, input->core.photoperiod);
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

    if (input->core.light == REAL(0.0)) {
        if (!trigger) {
            starch_night_start = input->core.starch;
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
    simulation_result_init(result, cfg.days, input->core.photoperiod, cfg.dt_hours);
    /* TODO: Validate allocation results and return an explicit error code instead of assuming success. */

    input->core.lambda_sb = lambda_sb_f(
        input->core.lambda_sb,
        input->core.nitrogen_soil_content,
        input->core.phosphorus_soil_content);

    const ModelInterface* model = simulation_plant_model_interface();
    ModelInterface runtime_model = *model;
    ModelBuffers buffers = {0};
    ModelStatus status = MODEL_STATUS_OK;

    PlantCtx rhs_ctx = {
        .base = input,
        .starch_night_start = input->core.starch
    };
    runtime_model.model_ctx = &rhs_ctx;

    Input tmp = *input;
    update_light_conditions(&tmp, REAL(0.0));
    rhs_ctx.was_light = (tmp.core.light != REAL(0.0));

    SimulationRuntime runtime;
    simulation_prepare_runtime(&runtime, &cfg, input, result, &rhs_ctx);

    status = simulation_run_model(&cfg, input, &runtime_model, &buffers, &runtime);
    if (status == MODEL_STATUS_OK)
        model_buffers_free(&buffers);

    if (status != MODEL_STATUS_OK && buffers.state)
        model_buffers_free(&buffers);

    simulation_runtime_cleanup(&runtime);
}

void simulation_result_free(SimulationResult* result)
{
    if (!result)
        return;
    free(result->sucrose);
    free(result->starch);
    free(result->ph);
    free(result->partition);
    result->sucrose = NULL;
    result->starch = NULL;
    result->ph = NULL;
    result->partition = NULL;
}

real_t compute_rgr(real_t FW_start, real_t FW_end, real_t duration_hours)
{
    return (RLOG(FW_end) - RLOG(FW_start)) / duration_hours;
}

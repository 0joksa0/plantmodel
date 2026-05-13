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

static void simulation_sync_buffers_from_input(ModelBuffers* buffers, const Input* input)
{
    simulation_pack_state(buffers->state, input);
}

static void simulation_sync_input_from_buffers(Input* input, const ModelBuffers* buffers)
{
    simulation_unpack_state(input, buffers->state);
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
    SimulationConfig cfg = config ? *config : simulation_default_config();
    if (cfg.days <= 0)
        cfg.days = 1;
    if (cfg.dt_hours <= REAL(0.0))
        cfg.dt_hours = DEFAULT_SIMULATION_DT;

    int days = cfg.days;
    result->days = days;
    result->photoperiod = input->core.photoperiod;
    result->total_steps = (int)(days * REAL(24.0) / cfg.dt_hours);
    result->filled_steps = 0;

    result->sucrose = calloc((size_t)result->total_steps, sizeof(real_t));
    result->starch = calloc((size_t)result->total_steps, sizeof(real_t));
    result->ph = calloc((size_t)result->total_steps, sizeof(real_t));
    result->partition = calloc((size_t)result->total_steps, sizeof(real_t));
    /* TODO: Validate allocation results and return an explicit error code instead of assuming success. */

    input->core.lambda_sb = lambda_sb_f(
        input->core.lambda_sb,
        input->core.nitrogen_soil_content,
        input->core.phosphorus_soil_content);

    const ModelInterface* model = simulation_plant_model_interface();
    ModelInterface runtime_model = *model;
    ModelBuffers buffers = {0};
    ModelStatus status;

    PlantCtx rhs_ctx = {
        .base = input,
        .starch_night_start = input->core.starch
    };
    runtime_model.model_ctx = &rhs_ctx;

    Input tmp = *input;
    update_light_conditions(&tmp, REAL(0.0));
    rhs_ctx.was_light = (tmp.core.light != REAL(0.0));

    char filename[64];
    snprintf(filename, sizeof(filename), "plant_%dh.csv", (int)input->core.photoperiod);

    CSVObserverCtx csv_ctx = {
        .f = cfg.write_csv ? fopen(filename, "w") : NULL,
        .plant_ctx = &rhs_ctx,
        .unpack_state = simulation_unpack_state,
        .compute_algebraic = simulation_compute_algebraic
    };
    /* TODO: Surface CSV open failures to caller/tests (currently silently continues without output file). */

    if (csv_ctx.f) {
        fprintf(csv_ctx.f,
            "t,light,photosynthesis,"
            "leaf_biomass,root_biomass,total_biomass,"
            "sucrose,starch,"
            "nitrogen,phosphorus,"
            "nitrogen_uptake,phosphorus_uptake,"
            "uptake_cost,transport_cost,"
            "starch_degradation_rate,night_efficiency\n");

        fflush(csv_ctx.f);
    }

    SimObsCtx sample_ctx = {
        .next_t = REAL(0.0),
        .dt_sample = cfg.dt_hours,
        .step = 0,
        .max_steps = result->total_steps,
        .rhs_ctx = &rhs_ctx,
        .result = result,
        .unpack_state = simulation_unpack_state,
        .compute_algebraic = simulation_compute_algebraic
    };
    NightStartObserverCtx night_ctx = {
        .plant_ctx = &rhs_ctx,
        .starch_index = X_STARCH
    };

    VectorObserver obs_list[4] = {
        nan_guard_observer,
        night_start_observer,
        sim_sampling_observer,
        csv_ctx.f ? csv_observer : NULL
    };
    void* ctx_list[4] = {
        NULL,
        &night_ctx,
        &sample_ctx,
        csv_ctx.f ? (void*)&csv_ctx : NULL
    };
    ObserverChain chain = {
        .obs = obs_list,
        .ctx = ctx_list,
        .count = 4
    };
    SolverStats stats;
    solver_stats_init(&stats);

    SolverPluginManager plugins;
    solver_plugin_manager_init(&plugins);

    SolverLogger logger;
    SolverLogger* logger_ptr = NULL;
    if (cfg.write_csv) {
        if (solver_logger_init(&logger, "steps.csv", "summary.csv") == 0) {
            logger_ptr = &logger;
        }
    }

    SolverDefaultRuntimeConfig runtime_cfg = {
        .observer = observer_chain_cb,
        .observer_ctx = &chain,
        .stats = &stats,
        .logger = logger_ptr
    };
    if (solver_default_runtime_init(&plugins, &runtime_cfg) != 0) {
        /* Fallback: run without runtime plugins if default runtime setup fails. */
        plugins.count = 0;
    }
    GuiOutput gui_output = {
        .sucrose = result->sucrose,
        .starch = result->starch,
        .photosynthesis = result->ph,
        .total_biomass = result->partition,
        .filled_steps = &result->filled_steps,
        .total_steps = result->total_steps,
        .dt_hours = cfg.dt_hours,
        .photoperiod = result->photoperiod,
        .days = result->days
    };
    GuiThreadArgs gui_args = {
        .output = &gui_output,
        .title = "Plant Simulation",
        .target_fps = 60,
        .summary_csv_path = "summary.csv"
    };

    pthread_t gui_thread;
    int gui_started = 0;
    if (cfg.gui_enabled) {
        gui_started = (pthread_create(&gui_thread, NULL, gui_main_thread, &gui_args) == 0);
    }
    /* TODO: Surface GUI thread startup failures so caller can react (headless fallback/logging). */

    SolverRunConfig run_cfg = {
        .type = cfg.solver_type,
        .t0 = REAL(0.0),
        .t_end = (real_t)days * REAL(24.0),
        .h_init = cfg.dt_hours,
        .tol = cfg.solver_tolerance,
        .plugins = &plugins
    };

    status = model_buffers_init(&buffers, runtime_model.descriptor);
    if (status != MODEL_STATUS_OK) {
        fprintf(stderr, "model_buffers_init failed: %s\n", model_status_string(status));
        if (csv_ctx.f)
            fclose(csv_ctx.f);
        if (logger_ptr) {
            solver_logger_close(&logger);
        }
        if (gui_started) {
            pthread_join(gui_thread, NULL);
        }
        return;
    }

    simulation_sync_buffers_from_input(&buffers, input);
    status = model_run(&runtime_model, &buffers, &run_cfg);
    if (status != MODEL_STATUS_OK) {
        fprintf(stderr, "model_run failed: %s\n", model_status_string(status));
        model_buffers_free(&buffers);
        if (csv_ctx.f)
            fclose(csv_ctx.f);
        if (logger_ptr) {
            solver_logger_close(&logger);
        }
        if (gui_started) {
            pthread_join(gui_thread, NULL);
        }
        return;
    }

    simulation_sync_input_from_buffers(input, &buffers);
    model_buffers_free(&buffers);

    if (csv_ctx.f)
        fclose(csv_ctx.f);

    if (logger_ptr) {
        solver_logger_close(&logger);
    }

    if (gui_started)
        pthread_join(gui_thread, NULL);
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

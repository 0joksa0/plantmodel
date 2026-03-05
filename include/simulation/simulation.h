#ifndef SIMULATION_H
#define SIMULATION_H

#include "model/input.h"
#include <solver.h>

typedef struct SimulationResult {
    real_t* sucrose;
    real_t* starch;
    real_t* ph;
    real_t* partition;
    int filled_steps;
    int total_steps;
    real_t photoperiod;
    int days;
} SimulationResult;

typedef struct SimulationConfig {
    int days;
    real_t dt_hours;
    real_t solver_tolerance;
    SolverType solver_type;
    int gui_enabled;
    int write_csv;
} SimulationConfig;

SimulationConfig simulation_default_config(void);

/* TODO: Document ownership/lifetime of `SimulationResult` buffers and return status from `simulate_days`. */
void simulate_days(const SimulationConfig* config, Input* input, SimulationResult* result);

void simulation_result_free(SimulationResult* result);

real_t compute_rgr(real_t FW_start, real_t FW_end, real_t duration_hours);

#endif

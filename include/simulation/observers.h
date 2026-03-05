#ifndef SIMULATION_OBSERVERS_H
#define SIMULATION_OBSERVERS_H

#include "model/input.h"
#include "simulation/simulation.h"
#include <stddef.h>
#include <solver.h>
#include <stdio.h>

typedef struct {
    Input* base;
    real_t starch_night_start;
    int was_light;
} PlantCtx;

typedef struct {
    real_t next_t;
    real_t dt_sample;
    int step;
    int max_steps;
    PlantCtx* rhs_ctx;
    SimulationResult* result;
    void (*unpack_state)(Input* in, const real_t* x);
    void (*compute_algebraic)(Input* s, real_t hour_of_day, real_t starch_night_start);
} SimObsCtx;

typedef struct {
    FILE* f;
    PlantCtx* plant_ctx;
    void (*unpack_state)(Input* in, const real_t* x);
    void (*compute_algebraic)(Input* s, real_t hour_of_day, real_t starch_night_start);
} CSVObserverCtx;

typedef struct {
    PlantCtx* plant_ctx;
    size_t starch_index;
} NightStartObserverCtx;

typedef struct {
    VectorObserver* obs;
    void** ctx;
    size_t count;
} ObserverChain;

void observer_chain_cb(real_t t, const real_t* x, size_t n, void* vctx);
void sim_sampling_observer(real_t t, const real_t* x, size_t n, void* vctx);
void csv_observer(real_t t, const real_t* x, size_t n, void* ctx);
void nan_guard_observer(real_t t, const real_t* x, size_t n, void* ctx);
void night_start_observer(real_t t, const real_t* x, size_t n, void* vctx);

#endif

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

/*
 * Divergence record for nan_guard_observer.
 *
 * The guard used to call abort() on the first non-finite state. It is called
 * from batch runners (parameter sweeps, cross-validation) that must classify a
 * divergent parameter combination and carry on, so it now records the first
 * occurrence and warns once instead. Callers that need to know should test
 * `tripped` after simulate_days, or check their own outputs with pm_isfinite
 * (include/util/fpclass.h) — not isfinite, which -ffast-math folds away.
 */
typedef struct {
    int    tripped;      /* non-zero once a non-finite state has been seen */
    double t_first;      /* simulation time [h] of the first occurrence */
    size_t index_first;  /* state index (PlantStateIndex) that went non-finite */
    long   n_events;     /* total non-finite samples observed */
    int    warned;       /* stderr warning already emitted */
} NanGuardCtx;

/*
 * Admissible-domain record for domain_guard_observer.
 *
 * Distinct from NanGuardCtx, and needed for a different reason.  The nutrient
 * pools N and P are contents: the model is defined on N >= 0, P >= 0 and its
 * saturating ratios are only monotone there.  Nothing in the ODEs enforces
 * that, and under a high carbon input the uptake terms cannot keep up with
 * demand, so the pools cross zero and the trajectory continues into a region
 * where the equations no longer mean what they were written to mean.  The
 * states stay finite throughout, so the NaN guard is silent and the run looks
 * healthy to every existing check.
 *
 * Recording the crossing lets a caller tell the two failure modes apart:
 * a run that merely became stiff and hit the step budget (integrable, the
 * trajectory is physical, the budget was the binding constraint) from a run
 * whose state left the domain the model is valid on (the numbers are finite
 * but not meaningful).  Before this record existed both were reported as
 * `diverged`, which is why the divergent region was treated as infeasible
 * rather than as inadmissible-for-a-physical-reason.
 */
typedef struct {
    int    tripped;      /* non-zero once a pool has gone negative */
    double t_first;      /* simulation time [h] of the first crossing */
    size_t index_first;  /* state index (PlantStateIndex) that crossed first */
    double min_value;    /* most negative value reached by any guarded pool */
    size_t index_min;    /* state index attaining min_value */
    long   n_events;     /* sampled states with at least one negative pool */
} DomainGuardCtx;

typedef struct {
    VectorObserver* obs;
    void** ctx;
    size_t count;
} ObserverChain;

void observer_chain_cb(real_t t, const real_t* x, size_t n, void* vctx);
void sim_sampling_observer(real_t t, const real_t* x, size_t n, void* vctx);
void csv_observer(real_t t, const real_t* x, size_t n, void* ctx);
void nan_guard_observer(real_t t, const real_t* x, size_t n, void* ctx);
void domain_guard_observer(real_t t, const real_t* x, size_t n, void* ctx);
void night_start_observer(real_t t, const real_t* x, size_t n, void* vctx);

#endif

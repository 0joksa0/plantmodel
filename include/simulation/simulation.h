#ifndef SIMULATION_H
#define SIMULATION_H

#include "model/input.h"
#include <solver.h>

typedef struct SimulationResult {
    real_t* sucrose;          /* cytosolic sucrose [µmol/gFW] */
    real_t* starch;
    real_t* ph;
    /*
     * Total biomass [gFW], NOT the starch/sucrose partition coefficient.
     *
     * This buffer was called `partition` and every consumer in tools/ reads it
     * as biomass (validation_runner, horizon_study), which is what
     * sim_sampling_observer writes into it and what the GUI binds it to.  The
     * name invited the opposite reading and produced a reported partition
     * coefficient of 152 at 12 h -- that number was the plant's total fresh
     * weight on day 30.  The partition coefficient is state X_STARCH_PART and
     * is not sampled into SimulationResult at all; read it from the solver
     * state vector (see probe/kstrue.c) or from the StarchPartitionCoeff
     * column of the export CSV.
     */
    real_t* total_biomass;
    real_t* sucrose_vacuole;  /* vacuolar sucrose [µmol/gFW] */
    real_t* sucrose_total;    /* cytosol + vacuole — compared against experiments */
    int filled_steps;
    int total_steps;
    real_t photoperiod;
    int days;
    /*
     * Divergence report, filled by nan_guard_observer during the run.
     * `diverged` non-zero means at least one state became non-finite, so the
     * end-of-run Input and every sampled buffer are meaningless for this
     * parameter combination — the caller must classify the run as invalid
     * rather than reading a growth rate out of it. Batch runners should test
     * this field (see tools/multi_2d_sweep.c) instead of testing outputs with
     * isfinite, which -ffast-math folds to a constant.
     */
    int diverged;
    double diverged_t;      /* simulation time [h] of first non-finite state */
    int diverged_index;     /* PlantStateIndex that first went non-finite; -1 if none */

    /*
     * Why the run failed, when it failed.  `diverged` alone conflates three
     * different situations, and treating them alike is what made a whole
     * region of parameter space look infeasible:
     *
     *   nonfinite   the state blew up; nothing here is usable.
     *   stalled     the step budget ran out while every state stayed finite.
     *               The trajectory up to `diverged_t` is a valid solution; the
     *               run is stiff, not wrong.
     *   left_domain a nutrient pool went negative.  The model's saturating
     *               ratios are only defined for N >= 0, P >= 0, so beyond
     *               `domain_exit_t` the equations no longer describe the
     *               intended biology even though the numbers stay finite.
     *
     * The last two are independent: a run can leave the domain and still
     * integrate to t_end, and that case is now reported rather than hidden.
     */
    int stalled;            /* step budget exhausted with all states finite */
    int left_domain;        /* a nutrient pool went negative at some point */
    double domain_exit_t;   /* simulation time [h] of the first crossing */
    int domain_exit_index;  /* PlantStateIndex that crossed first; -1 if none */
    double domain_min;      /* most negative pool value reached */
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

#ifndef SIMULATION_DYNAMICS_H
#define SIMULATION_DYNAMICS_H

#include "model/input.h"
#include "model/state.h"
#include "simulation/observers.h"
#include <framework/model_runner.h>
#include <solver.h>

void simulation_pack_state(real_t* x, const Input* in);
void simulation_unpack_state(Input* in, const real_t* x);
void simulation_compute_algebraic(Input* s, real_t hour_of_day, real_t starch_night_start);
void simulation_plant_rhs(real_t t, const real_t* x, real_t* dx, void* user_ctx);
void simulation_plant_rhs_adapter(real_t t, const real_t* state, real_t* dst, size_t n, void* ctx);
const ModelInterface* simulation_plant_model_interface(void);

#endif

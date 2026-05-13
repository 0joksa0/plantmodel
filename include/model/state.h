#ifndef MODEL_STATE_H
#define MODEL_STATE_H

#include "model/input.h"
#include <solver.h>

typedef enum {
    X_STARCH_PART = 0,
    X_STARCH,
    X_SUCROSE,
    X_N_AFF,
    X_P_AFF,
    X_N,
    X_P,
    X_SUCROSE_ROOT_ALLOC,
    X_LEAF_BIOMASS,
    X_ROOT_BIOMASS,
    X_DIM
} PlantStateIndex;

typedef struct PlantState {
    real_t starch_partition_coeff;
    real_t starch;
    real_t sucrose;
    real_t nitrogen_affinity;
    real_t phosphorus_affinity;
    real_t nitrogen;
    real_t phosphorus;
    real_t sucrose_root_allocation;
    real_t leaf_biomass;
    real_t root_biomass;
} PlantState;

PlantState plant_state_from_input(const Input* input);
void plant_state_apply_to_input(Input* input, const PlantState* state);
void plant_state_pack(real_t* dst, const PlantState* state);
void plant_state_unpack(PlantState* dst, const real_t* src);

#endif

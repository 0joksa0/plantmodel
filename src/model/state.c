#include "model/state.h"

PlantState plant_state_from_input(const Input* input)
{
    PlantState state = {
        .starch_partition_coeff = input->carbohydrates.starch_partition_coeff,
        .starch = input->carbohydrates.starch,
        .sucrose = input->carbohydrates.sucrose,
        .nitrogen_affinity = input->nutrients.nitrogen_affinity,
        .phosphorus_affinity = input->nutrients.phosphorus_affinity,
        .nitrogen = input->nutrients.nitrogen,
        .phosphorus = input->nutrients.phosphorus,
        .sucrose_root_allocation = input->growth.sucrose_root_allocation,
        .leaf_biomass = input->growth.leaf_biomass,
        .root_biomass = input->growth.root_biomass
    };
    return state;
}

void plant_state_apply_to_input(Input* input, const PlantState* state)
{
    input->carbohydrates.starch_partition_coeff = state->starch_partition_coeff;
    input->carbohydrates.starch = state->starch;
    input->carbohydrates.sucrose = state->sucrose;
    input->nutrients.nitrogen_affinity = state->nitrogen_affinity;
    input->nutrients.phosphorus_affinity = state->phosphorus_affinity;
    input->nutrients.nitrogen = state->nitrogen;
    input->nutrients.phosphorus = state->phosphorus;
    input->growth.sucrose_root_allocation = state->sucrose_root_allocation;
    input->growth.leaf_biomass = state->leaf_biomass;
    input->growth.root_biomass = state->root_biomass;
    input->growth.total_biomass = input->growth.leaf_biomass + input->growth.root_biomass;
}

void plant_state_pack(real_t* dst, const PlantState* state)
{
    dst[X_STARCH_PART] = state->starch_partition_coeff;
    dst[X_STARCH] = state->starch;
    dst[X_SUCROSE] = state->sucrose;
    dst[X_N_AFF] = state->nitrogen_affinity;
    dst[X_P_AFF] = state->phosphorus_affinity;
    dst[X_N] = state->nitrogen;
    dst[X_P] = state->phosphorus;
    dst[X_SUCROSE_ROOT_ALLOC] = state->sucrose_root_allocation;
    dst[X_LEAF_BIOMASS] = state->leaf_biomass;
    dst[X_ROOT_BIOMASS] = state->root_biomass;
}

void plant_state_unpack(PlantState* dst, const real_t* src)
{
    dst->starch_partition_coeff = src[X_STARCH_PART];
    dst->starch = src[X_STARCH];
    dst->sucrose = src[X_SUCROSE];
    dst->nitrogen_affinity = src[X_N_AFF];
    dst->phosphorus_affinity = src[X_P_AFF];
    dst->nitrogen = src[X_N];
    dst->phosphorus = src[X_P];
    dst->sucrose_root_allocation = src[X_SUCROSE_ROOT_ALLOC];
    dst->leaf_biomass = src[X_LEAF_BIOMASS];
    dst->root_biomass = src[X_ROOT_BIOMASS];
}

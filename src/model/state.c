#include "model/state.h"

PlantState plant_state_from_input(const Input* input)
{
    PlantState state = {
        .starch_partition_coeff = input->core.starch_partition_coeff,
        .starch = input->core.starch,
        .sucrose = input->core.sucrose,
        .nitrogen_affinity = input->core.nitrogen_affinity,
        .phosphorus_affinity = input->core.phosphorus_affinity,
        .nitrogen = input->core.nitrogen,
        .phosphorus = input->core.phosphorus,
        .sucrose_root_allocation = input->core.sucrose_root_allocation,
        .leaf_biomass = input->core.leaf_biomass,
        .root_biomass = input->core.root_biomass
    };
    return state;
}

void plant_state_apply_to_input(Input* input, const PlantState* state)
{
    input->core.starch_partition_coeff = state->starch_partition_coeff;
    input->core.starch = state->starch;
    input->core.sucrose = state->sucrose;
    input->core.nitrogen_affinity = state->nitrogen_affinity;
    input->core.phosphorus_affinity = state->phosphorus_affinity;
    input->core.nitrogen = state->nitrogen;
    input->core.phosphorus = state->phosphorus;
    input->core.sucrose_root_allocation = state->sucrose_root_allocation;
    input->core.leaf_biomass = state->leaf_biomass;
    input->core.root_biomass = state->root_biomass;
    input->core.total_biomass = input->core.leaf_biomass + input->core.root_biomass;
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

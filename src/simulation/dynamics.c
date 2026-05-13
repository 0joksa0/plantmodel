#include "simulation/dynamics.h"

#include "model/algebraic.h"
#include "model/model.h"
#include "model/outputs.h"
#include "model/parameters.h"

#include <math.h>

static const ModelFieldDescriptor PLANT_STATE_FIELDS[] = {
    { .name = "starch_partition_coeff", .kind = MODEL_FIELD_STATE, .index = X_STARCH_PART, .unit = NULL, .description = "Starch partition coefficient" },
    { .name = "starch", .kind = MODEL_FIELD_STATE, .index = X_STARCH, .unit = NULL, .description = "Starch pool" },
    { .name = "sucrose", .kind = MODEL_FIELD_STATE, .index = X_SUCROSE, .unit = NULL, .description = "Sucrose pool" },
    { .name = "nitrogen_affinity", .kind = MODEL_FIELD_STATE, .index = X_N_AFF, .unit = NULL, .description = "Nitrogen affinity" },
    { .name = "phosphorus_affinity", .kind = MODEL_FIELD_STATE, .index = X_P_AFF, .unit = NULL, .description = "Phosphorus affinity" },
    { .name = "nitrogen", .kind = MODEL_FIELD_STATE, .index = X_N, .unit = NULL, .description = "Nitrogen content" },
    { .name = "phosphorus", .kind = MODEL_FIELD_STATE, .index = X_P, .unit = NULL, .description = "Phosphorus content" },
    { .name = "sucrose_root_allocation", .kind = MODEL_FIELD_STATE, .index = X_SUCROSE_ROOT_ALLOC, .unit = NULL, .description = "Sucrose root allocation" },
    { .name = "leaf_biomass", .kind = MODEL_FIELD_STATE, .index = X_LEAF_BIOMASS, .unit = NULL, .description = "Leaf biomass" },
    { .name = "root_biomass", .kind = MODEL_FIELD_STATE, .index = X_ROOT_BIOMASS, .unit = NULL, .description = "Root biomass" }
};

static const ModelDescriptor PLANT_MODEL_DESCRIPTOR = {
    .name = "plantmodel",
    .state_count = X_DIM,
    .parameter_count = 0,
    .input_count = 0,
    .output_count = 0,
    .state_fields = PLANT_STATE_FIELDS,
    .parameter_fields = NULL,
    .input_fields = NULL,
    .output_fields = NULL
};

static void simulation_plant_model_rhs(
    real_t t,
    const real_t* state,
    const real_t* parameters,
    const real_t* inputs,
    real_t* dstate,
    void* model_ctx)
{
    (void)parameters;
    (void)inputs;
    simulation_plant_rhs(t, state, dstate, model_ctx);
}

static const ModelInterface PLANT_MODEL_INTERFACE = {
    .descriptor = &PLANT_MODEL_DESCRIPTOR,
    .callbacks = {
        .rhs = simulation_plant_model_rhs,
        .outputs = NULL
    },
    .model_ctx = NULL
};

void simulation_pack_state(real_t* x, const Input* in)
{
    PlantState state = plant_state_from_input(in);
    plant_state_pack(x, &state);
}

void simulation_unpack_state(Input* in, const real_t* x)
{
    PlantState state = {0};
    plant_state_unpack(&state, x);
    plant_state_apply_to_input(in, &state);
}

void simulation_compute_algebraic(
    Input* s,
    real_t hour_of_day,
    real_t starch_night_start)
{
    PlantState state = plant_state_from_input(s);
    PlantParameters parameters = plant_parameters_from_input(s);
    PlantEnvironment environment = {
        .hour_of_day = hour_of_day
    };
    PlantOutputs outputs = {0};

    plant_compute_outputs(&state, &parameters, &environment, starch_night_start, &outputs);
    plant_outputs_apply_to_input(s, &outputs);
}

void simulation_plant_rhs(
    real_t t,
    const real_t* x,
    real_t* dx,
    void* user_ctx)
{
    PlantCtx* ctx = (PlantCtx*)user_ctx;

    Input s = *(ctx->base);

    s.core.starch_partition_coeff = x[X_STARCH_PART];
    s.core.starch = x[X_STARCH];
    s.core.sucrose = x[X_SUCROSE];
    s.core.nitrogen_affinity = x[X_N_AFF];
    s.core.phosphorus_affinity = x[X_P_AFF];
    s.core.nitrogen = x[X_N];
    s.core.phosphorus = x[X_P];
    s.core.sucrose_root_allocation = x[X_SUCROSE_ROOT_ALLOC];
    s.core.leaf_biomass = x[X_LEAF_BIOMASS];
    s.core.root_biomass = x[X_ROOT_BIOMASS];
    s.core.total_biomass = s.core.leaf_biomass + s.core.root_biomass;

    real_t hour_of_day = fmod(t, REAL(24.0));

    simulation_compute_algebraic(&s, hour_of_day, ctx->starch_night_start);

    dx[X_STARCH_PART] = starch_sucrose_partition(
        s.core.light,
        s.core.starch_partition_coeff,
        s.core.min_sucrose,
        s.core.max_sucrose,
        s.core.sucrose,
        s.core.lambda_sdr,
        s.core.lambda_sdi,
        s.core.lambda_sni);

    dx[X_STARCH] = starch_production(
        s.core.photosynthesis,
        s.core.starch_partition_coeff,
        s.core.starch_degradation_rate);

    dx[X_SUCROSE] = sucrose_production(
        s.core.starch_partition_coeff,
        s.core.photosynthesis,
        s.core.starch_degradation_rate,
        s.core.uptake_cost,
        s.core.transport_cost,
        s.core.respiration_frequency,
        s.core.sucrose,
        s.core.sucrose_loading_frequency,
        s.core.night_efficiency_starch);

    dx[X_N_AFF] = nitrogen_affinity(
        s.core.nitrogen_affinity,
        s.core.nitrogen_uptake,
        s.core.nitrogen_cost,
        s.core.max_nitrogen,
        s.core.nitrogen,
        s.core.photosynthesis,
        s.core.max_photosynthetic_rate,
        s.core.lambda_k,
        s.core.optimal_stoichiometric_ratio,
        s.core.phosphorus,
        s.core.min_nitrogen,
        s.core.nitrogen_uptake_sucrose_consumption,
        s.core.starch_partition_coeff,
        s.core.starch_degradation_rate);

    dx[X_P_AFF] = phosphorus_affinity(
        s.core.phosphorus_affinity,
        s.core.phosphorus_uptake,
        s.core.phosphorus_cost,
        s.core.max_phosphorus,
        s.core.phosphorus,
        s.core.photosynthesis,
        s.core.max_photosynthetic_rate,
        s.core.lambda_k,
        s.core.optimal_stoichiometric_ratio,
        s.core.nitrogen,
        s.core.min_phosphorus,
        s.core.phosphorus_uptake_sucrose_consumption,
        s.core.starch_partition_coeff,
        s.core.starch_degradation_rate);

    dx[X_N] = nitrogen_content(
        s.core.nitrogen_uptake,
        s.core.respiration_frequency,
        s.core.sucrose,
        s.core.starch,
        s.core.sucrose_loading_frequency,
        s.core.night_efficiency_starch,
        s.core.assimilation_cost_nitrogen,
        s.core.leaf_biomass,
        s.core.total_biomass,
        s.core.photosynthesis,
        s.core.max_photosynthetic_rate,
        s.core.nutrient_conversion_parameter,
        s.core.min_nitrogen_photosynthesis);

    dx[X_P] = phosphorus_content(
        s.core.phosphorus_uptake,
        s.core.respiration_frequency,
        s.core.sucrose,
        s.core.starch,
        s.core.sucrose_loading_frequency,
        s.core.night_efficiency_starch,
        s.core.assimilation_cost_phosphorus,
        s.core.leaf_biomass,
        s.core.total_biomass,
        s.core.photosynthesis,
        s.core.max_photosynthetic_rate,
        s.core.nutrient_conversion_parameter,
        s.core.min_phosphorus_photosynthesis);

    dx[X_SUCROSE_ROOT_ALLOC] = sucrose_root_allocation(
        s.core.sucrose_root_allocation,
        s.core.nitrogen_affinity,
        s.core.phosphorus_affinity,
        s.core.stoichiometric_signal,
        s.core.sucrose,
        s.core.min_sucrose,
        s.core.nitrogen,
        s.core.min_nitrogen,
        s.core.phosphorus,
        s.core.min_phosphorus);

    dx[X_LEAF_BIOMASS] = leaf_growth(
        s.core.lambda_sb,
        s.core.sucrose_root_allocation,
        s.core.sucrose_loading_frequency,
        s.core.night_efficiency_starch,
        s.core.leaf_biomass,
        s.core.leaf_deathrate,
        s.core.leaf_competitive_rate);

    dx[X_ROOT_BIOMASS] = root_growth(
        s.core.lambda_sb,
        s.core.sucrose_root_allocation,
        s.core.sucrose_loading_frequency,
        s.core.night_efficiency_starch,
        s.core.leaf_biomass,
        s.core.root_biomass,
        s.core.root_deathrate,
        s.core.root_competitive_rate);
}

void simulation_plant_rhs_adapter(
    real_t t,
    const real_t* state,
    real_t* dst,
    size_t n,
    void* ctx)
{
    (void)n;
    simulation_plant_rhs(t, state, dst, ctx);
}

const ModelInterface* simulation_plant_model_interface(void)
{
    return &PLANT_MODEL_INTERFACE;
}

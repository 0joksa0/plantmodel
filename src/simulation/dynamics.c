#include "simulation/dynamics.h"

#include "model/algebraic.h"
#include "model/model.h"
#include "model/outputs.h"
#include "model/parameters.h"

#include <math.h>

static const ModelFieldDescriptor PLANT_STATE_FIELDS[] = {
    { .name = "starch_partition_coeff", .kind = MODEL_FIELD_STATE, .index = X_STARCH_PART,        .unit = NULL, .description = "Starch partition coefficient" },
    { .name = "starch",                 .kind = MODEL_FIELD_STATE, .index = X_STARCH,             .unit = NULL, .description = "Starch pool" },
    { .name = "sucrose",                .kind = MODEL_FIELD_STATE, .index = X_SUCROSE,            .unit = NULL, .description = "Cytosolic sucrose pool [umol/gFW]" },
    { .name = "nitrogen_affinity",      .kind = MODEL_FIELD_STATE, .index = X_N_AFF,              .unit = NULL, .description = "Nitrogen affinity" },
    { .name = "phosphorus_affinity",    .kind = MODEL_FIELD_STATE, .index = X_P_AFF,              .unit = NULL, .description = "Phosphorus affinity" },
    { .name = "nitrogen",               .kind = MODEL_FIELD_STATE, .index = X_N,                  .unit = NULL, .description = "Nitrogen content" },
    { .name = "phosphorus",             .kind = MODEL_FIELD_STATE, .index = X_P,                  .unit = NULL, .description = "Phosphorus content" },
    { .name = "sucrose_root_allocation",.kind = MODEL_FIELD_STATE, .index = X_SUCROSE_ROOT_ALLOC, .unit = NULL, .description = "Sucrose root allocation" },
    { .name = "leaf_biomass",           .kind = MODEL_FIELD_STATE, .index = X_LEAF_BIOMASS,       .unit = NULL, .description = "Leaf biomass" },
    { .name = "root_biomass",           .kind = MODEL_FIELD_STATE, .index = X_ROOT_BIOMASS,       .unit = NULL, .description = "Root biomass" },
    { .name = "sucrose_vacuole",        .kind = MODEL_FIELD_STATE, .index = X_SUCROSE_VAC,        .unit = NULL, .description = "Vacuolar sucrose pool [umol/gFW]" }
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

    /* Pool contents enter fifteen saturating ratios of the form
     * pool/(pool + K).  Each has a pole at pool = -K, and the poles sit at
     * small negative values: max_nitrogen, min_nitrogen, max_sucrose and so on
     * are all O(1) to O(100) in these units.  A Runge-Kutta trial stage that
     * takes a depleted pool slightly negative therefore does not merely give a
     * slightly wrong rate -- it can land arbitrarily close to a pole and return
     * a derivative many orders of magnitude too large.
     *
     * This is not hypothetical: with sucrose exhausted the nitrogen pool drains
     * (uptake is gated on sucrose and goes to zero, while the respiratory debit
     * keeps drawing on starch), crosses zero, meets the pole at N = -max_N in
     * nitrogen_affinity, and the state reaches 1e203 within one step.
     *
     * Clamping the pool contents to their non-negative part as they enter the
     * constitutive laws is the identity on the whole physically meaningful
     * region -- a pool content is by definition non-negative -- so no valid
     * trajectory is altered.  The state vector itself is left untouched, so
     * conservation is unaffected and the solver still sees the true state. */
    s.carbohydrates.starch_partition_coeff = x[X_STARCH_PART];
    s.carbohydrates.starch = RMAX(x[X_STARCH], REAL(0.0));
    s.carbohydrates.sucrose = RMAX(x[X_SUCROSE], REAL(0.0));
    s.nutrients.nitrogen_affinity = x[X_N_AFF];
    s.nutrients.phosphorus_affinity = x[X_P_AFF];
    s.nutrients.nitrogen = RMAX(x[X_N], REAL(0.0));
    s.nutrients.phosphorus = RMAX(x[X_P], REAL(0.0));
    s.growth.sucrose_root_allocation = x[X_SUCROSE_ROOT_ALLOC];
    s.growth.leaf_biomass = RMAX(x[X_LEAF_BIOMASS], REAL(0.0));
    s.growth.root_biomass = RMAX(x[X_ROOT_BIOMASS], REAL(0.0));
    s.growth.total_biomass = s.growth.leaf_biomass + s.growth.root_biomass;
    s.carbohydrates.sucrose_vacuole = RMAX(x[X_SUCROSE_VAC], REAL(0.0));

    real_t hour_of_day = fmod(t, REAL(24.0));

    simulation_compute_algebraic(&s, hour_of_day, ctx->starch_night_start);

    dx[X_STARCH_PART] = starch_sucrose_partition(
        s.photo.light,
        s.carbohydrates.starch_partition_coeff,
        s.carbohydrates.min_sucrose,
        s.carbohydrates.max_sucrose,
        s.carbohydrates.sucrose,
        s.carbohydrates.lambda_sdr,
        s.carbohydrates.lambda_sdi,
        s.carbohydrates.lambda_sni);

    dx[X_STARCH] = starch_production(
        s.photo.photosynthesis,
        s.carbohydrates.starch_partition_coeff,
        s.carbohydrates.starch_degradation_rate);

    /* Cytosolic sucrose: mass-conserving coupling with vacuolar pool.
     * At equilibrium k_in*s_cyt = k_out*s_vac, so these terms cancel and
     * the cytosolic steady-state is identical to the single-compartment model.
     * Ref: Lohaus et al. (2007) Planta 227:143-156; Kaiser & Heber (1984) Planta 161:562 */
    dx[X_SUCROSE] = sucrose_production(
        s.carbohydrates.starch_partition_coeff,
        s.photo.photosynthesis,
        s.carbohydrates.starch_degradation_rate,
        s.carbohydrates.uptake_cost,
        s.carbohydrates.transport_cost,
        s.carbohydrates.respiration_frequency,
        s.carbohydrates.sucrose,
        s.carbohydrates.sucrose_loading_frequency,
        s.carbohydrates.night_efficiency_starch,
        s.photo.light)
        - s.carbohydrates.k_vacuole_in  * s.carbohydrates.sucrose
        + s.carbohydrates.k_vacuole_out * s.carbohydrates.sucrose_vacuole;
    /* The construction-cost debit is subtracted after the two biomass
     * derivatives are assigned, at the end of this function. */

    dx[X_SUCROSE_VAC] = s.carbohydrates.k_vacuole_in  * s.carbohydrates.sucrose
                      - s.carbohydrates.k_vacuole_out * s.carbohydrates.sucrose_vacuole;

    dx[X_N_AFF] = nitrogen_affinity(
        s.nutrients.nitrogen_affinity,
        s.nutrients.nitrogen_uptake,
        s.nutrients.nitrogen_cost,
        s.nutrients.max_nitrogen,
        s.nutrients.nitrogen,
        s.photo.photosynthesis,
        s.photo.max_photosynthetic_rate,
        s.nutrients.lambda_k,
        s.nutrients.optimal_stoichiometric_ratio,
        s.nutrients.phosphorus,
        s.nutrients.min_nitrogen,
        s.nutrients.nitrogen_uptake_sucrose_consumption,
        s.carbohydrates.starch_partition_coeff,
        s.carbohydrates.starch_degradation_rate);

    dx[X_P_AFF] = phosphorus_affinity(
        s.nutrients.phosphorus_affinity,
        s.nutrients.phosphorus_uptake,
        s.nutrients.phosphorus_cost,
        s.nutrients.max_phosphorus,
        s.nutrients.phosphorus,
        s.photo.photosynthesis,
        s.photo.max_photosynthetic_rate,
        s.nutrients.lambda_k,
        s.nutrients.optimal_stoichiometric_ratio,
        s.nutrients.nitrogen,
        s.nutrients.min_phosphorus,
        s.nutrients.phosphorus_uptake_sucrose_consumption,
        s.carbohydrates.starch_partition_coeff,
        s.carbohydrates.starch_degradation_rate);

    dx[X_N] = nitrogen_content(
        s.nutrients.nitrogen_uptake,
        s.carbohydrates.respiration_frequency,
        s.carbohydrates.sucrose,
        s.carbohydrates.starch,
        s.carbohydrates.sucrose_loading_frequency,
        s.carbohydrates.night_efficiency_starch,
        s.nutrients.assimilation_cost_nitrogen,
        s.growth.leaf_biomass,
        s.growth.total_biomass,
        s.photo.photosynthesis,
        s.photo.max_photosynthetic_rate,
        s.nutrients.nutrient_conversion_parameter,
        s.nutrients.min_nitrogen_photosynthesis);

    dx[X_P] = phosphorus_content(
        s.nutrients.phosphorus_uptake,
        s.carbohydrates.respiration_frequency,
        s.carbohydrates.sucrose,
        s.carbohydrates.starch,
        s.carbohydrates.sucrose_loading_frequency,
        s.carbohydrates.night_efficiency_starch,
        s.nutrients.assimilation_cost_phosphorus,
        s.growth.leaf_biomass,
        s.growth.total_biomass,
        s.photo.photosynthesis,
        s.photo.max_photosynthetic_rate,
        s.nutrients.nutrient_conversion_parameter,
        s.nutrients.min_phosphorus_photosynthesis);

    dx[X_SUCROSE_ROOT_ALLOC] = sucrose_root_allocation(
        s.growth.sucrose_root_allocation,
        s.nutrients.nitrogen_affinity,
        s.nutrients.phosphorus_affinity,
        s.nutrients.stoichiometric_signal,
        s.carbohydrates.sucrose,
        s.carbohydrates.min_sucrose,
        s.nutrients.nitrogen,
        s.nutrients.min_nitrogen,
        s.nutrients.phosphorus,
        s.nutrients.min_phosphorus);

    dx[X_LEAF_BIOMASS] = leaf_growth(
        s.growth.lambda_sb,
        s.growth.sucrose_root_allocation,
        s.carbohydrates.sucrose_loading_frequency,
        s.carbohydrates.night_efficiency_starch,
        s.growth.leaf_biomass,
        s.growth.leaf_deathrate,
        s.growth.leaf_competitive_rate);

    dx[X_ROOT_BIOMASS] = root_growth(
        s.growth.lambda_sb,
        s.growth.sucrose_root_allocation,
        s.carbohydrates.sucrose_loading_frequency,
        s.carbohydrates.night_efficiency_starch,
        s.growth.leaf_biomass,
        s.growth.root_biomass,
        s.growth.root_deathrate,
        s.growth.root_competitive_rate);

    /* Construction cost of the biomass credited just above.  Growth is driven
     * by the phloem-load flux but, at HEAD, was never charged against it, so
     * the carbon budget closed identically while biomass grew freely.  The
     * stoichiometry is measured (growth respiration, Penning de Vries et al.),
     * not fitted; see growth_cost() in src/model/carbohydrates.c. */
    dx[X_SUCROSE] -= growth_cost(s.growth.growth_stoichiometry,
                                 dx[X_LEAF_BIOMASS], dx[X_ROOT_BIOMASS]);
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

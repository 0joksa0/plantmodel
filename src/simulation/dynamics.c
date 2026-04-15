#include "simulation/dynamics.h"

#include "model/model.h"

#include <math.h>

void simulation_pack_state(real_t* x, const Input* in)
{
    x[X_STARCH_PART] = in->carbohydrates.starch_partition_coeff;
    x[X_STARCH] = in->carbohydrates.starch;
    x[X_SUCROSE] = in->carbohydrates.sucrose;
    x[X_N_AFF] = in->nutrients.nitrogen_affinity;
    x[X_P_AFF] = in->nutrients.phosphorus_affinity;
    x[X_N] = in->nutrients.nitrogen;
    x[X_P] = in->nutrients.phosphorus;
    x[X_SUCROSE_ROOT_ALLOC] = in->growth.sucrose_root_allocation;
    x[X_LEAF_BIOMASS] = in->growth.leaf_biomass;
    x[X_ROOT_BIOMASS] = in->growth.root_biomass;
}

void simulation_unpack_state(Input* in, const real_t* x)
{
    in->carbohydrates.starch_partition_coeff = x[X_STARCH_PART];
    in->carbohydrates.starch = x[X_STARCH];
    in->carbohydrates.sucrose = x[X_SUCROSE];
    in->nutrients.nitrogen_affinity = x[X_N_AFF];
    in->nutrients.phosphorus_affinity = x[X_P_AFF];
    in->nutrients.nitrogen = x[X_N];
    in->nutrients.phosphorus = x[X_P];
    in->growth.sucrose_root_allocation = x[X_SUCROSE_ROOT_ALLOC];
    in->growth.leaf_biomass = x[X_LEAF_BIOMASS];
    in->growth.root_biomass = x[X_ROOT_BIOMASS];
    in->growth.total_biomass = in->growth.leaf_biomass + in->growth.root_biomass;
}

void simulation_compute_algebraic(
    Input* s,
    real_t hour_of_day,
    real_t starch_night_start)
{
    update_light_conditions(s, hour_of_day);

    s->nutrients.min_nitrogen = min_nitrogen(
        s->carbohydrates.respiration_frequency, s->carbohydrates.max_sucrose, s->carbohydrates.max_starch,
        s->carbohydrates.sucrose_loading_frequency, s->nutrients.assimilation_cost_nitrogen,
        s->nutrients.min_nitrogen_photosynthesis, s->nutrients.nutrient_conversion_parameter, s->photo.photoperiod);

    s->nutrients.min_phosphorus = min_phosphorus(
        s->carbohydrates.respiration_frequency, s->carbohydrates.max_sucrose, s->carbohydrates.max_starch,
        s->carbohydrates.sucrose_loading_frequency, s->nutrients.assimilation_cost_phosphorus,
        s->nutrients.min_phosphorus_photosynthesis, s->nutrients.nutrient_conversion_parameter, s->photo.photoperiod);

    s->carbohydrates.max_starch = max_starch(s->carbohydrates.max_starch_degradation_rate, s->photo.photoperiod);
    s->nutrients.max_nitrogen = max_nitrogen(s->nutrients.min_nitrogen);
    s->nutrients.max_phosphorus = max_phosphorus(s->nutrients.max_nitrogen, s->nutrients.optimal_stoichiometric_ratio);

    s->photo.nitrogen_saturation = nitrogen_saturation(s->nutrients.nitrogen, s->nutrients.min_nitrogen_photosynthesis);
    s->photo.phosphorus_saturation = phosphorus_saturation(
        s->nutrients.min_nitrogen_photosynthesis, s->nutrients.optimal_stoichiometric_ratio, s->nutrients.phosphorus);

    s->photo.limitation_of_photosynthetic_rate = limitation_of_photosynthetic_rate(s->carbohydrates.starch, s->photo.feedback_on_photosynthesis, s->carbohydrates.max_starch);

    s->photo.photosynthesis = farquhar_photosynthesis(s);
    // s->photo.photosynthesis = photosynthesis(s->photo.light, s->photo.limitation_of_photosynthetic_rate, s->photo.max_photosynthetic_rate, s->photo.nitrogen_saturation, s->photo.phosphorus_saturation, s->growth.leaf_biomass, s->photo.min_leaf_biomass);

    s->carbohydrates.night_efficiency_starch = night_efficiency_starch(
        s->carbohydrates.sucrose, s->carbohydrates.max_sucrose, s->carbohydrates.lambda_g, s->photo.light, s->carbohydrates.starch_partition_coeff);

    s->carbohydrates.starch_degradation_rate = starch_degradation(
        s->carbohydrates.max_starch_degradation_rate,
        s->carbohydrates.max_sucrose,
        starch_night_start,
        s->carbohydrates.min_starch,
        s->carbohydrates.sucrose,
        s->photo.light,
        hour_of_day,
        s->photo.photoperiod);

    s->carbohydrates.uptake_cost = uptake_cost(
        s->nutrients.nitrogen_uptake_sucrose_consumption, s->nutrients.nitrogen_uptake,
        s->nutrients.phosphorus_uptake_sucrose_consumption, s->nutrients.phosphorus_uptake);

    s->carbohydrates.transport_cost = transport_cost(
        s->carbohydrates.sucrose_consumption_transport, s->carbohydrates.respiration_frequency,
        s->carbohydrates.sucrose, s->carbohydrates.sucrose_loading_frequency, s->carbohydrates.night_efficiency_starch);

    s->nutrients.nitrogen_cost = nitrogen_cost(
        s->carbohydrates.respiration_frequency, s->carbohydrates.sucrose, s->carbohydrates.starch,
        s->carbohydrates.sucrose_loading_frequency, s->carbohydrates.night_efficiency_starch,
        s->nutrients.assimilation_cost_nitrogen,
        s->growth.leaf_biomass, s->growth.total_biomass,
        s->photo.photosynthesis, s->photo.max_photosynthetic_rate,
        s->nutrients.min_nitrogen_photosynthesis);

    s->nutrients.phosphorus_cost = phosphorus_cost(
        s->carbohydrates.respiration_frequency, s->carbohydrates.sucrose, s->carbohydrates.starch,
        s->carbohydrates.sucrose_loading_frequency, s->carbohydrates.night_efficiency_starch,
        s->nutrients.assimilation_cost_phosphorus,
        s->growth.leaf_biomass, s->growth.total_biomass,
        s->photo.photosynthesis, s->photo.max_photosynthetic_rate,
        s->nutrients.min_phosphorus_photosynthesis);

    s->nutrients.nitrogen_nutrient_uptake = nitrogen_nutrient_uptake(
        s->nutrients.max_nitrogen_uptake, s->nutrients.nitrogen_soil_content, s->nutrients.Michaelis_Menten_constant_nitrogen);

    s->nutrients.phosphorus_nutrient_uptake = phosphorus_nutrient_uptake(
        s->nutrients.max_phosphorus_uptake, s->nutrients.phosphorus_soil_content, s->nutrients.Michaelis_Menten_constant_phosphorus);

    s->nutrients.pot_nitrogen_uptake = pot_nitrogen_uptake(
        s->nutrients.nitrogen_nutrient_uptake, s->nutrients.nitrogen_affinity, s->growth.root_biomass, s->growth.total_biomass);

    s->nutrients.pot_phosphorus_uptake = pot_phosphorus_uptake(
        s->nutrients.phosphorus_nutrient_uptake, s->nutrients.phosphorus_affinity, s->growth.root_biomass, s->growth.total_biomass);

    s->nutrients.nitrogen_uptake = nitrogen_uptake(
        s->nutrients.pot_nitrogen_uptake, s->carbohydrates.sucrose, s->nutrients.nitrogen_uptake_sucrose_consumption);

    s->nutrients.phosphorus_uptake = phosphorus_uptake(
        s->nutrients.pot_phosphorus_uptake, s->carbohydrates.sucrose, s->nutrients.phosphorus_uptake_sucrose_consumption);

    s->nutrients.stoichiometric_signal = stoichiometric_signal(
        s->nutrients.optimal_stoichiometric_ratio, s->nutrients.nitrogen, s->nutrients.phosphorus);
}

void simulation_plant_rhs(
    real_t t,
    const real_t* x,
    real_t* dx,
    void* user_ctx)
{
    PlantCtx* ctx = (PlantCtx*)user_ctx;

    Input s = *(ctx->base);

    s.carbohydrates.starch_partition_coeff = x[X_STARCH_PART];
    s.carbohydrates.starch = x[X_STARCH];
    s.carbohydrates.sucrose = x[X_SUCROSE];
    s.nutrients.nitrogen_affinity = x[X_N_AFF];
    s.nutrients.phosphorus_affinity = x[X_P_AFF];
    s.nutrients.nitrogen = x[X_N];
    s.nutrients.phosphorus = x[X_P];
    s.growth.sucrose_root_allocation = x[X_SUCROSE_ROOT_ALLOC];
    s.growth.leaf_biomass = x[X_LEAF_BIOMASS];
    s.growth.root_biomass = x[X_ROOT_BIOMASS];
    s.growth.total_biomass = s.growth.leaf_biomass + s.growth.root_biomass;

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

    dx[X_SUCROSE] = sucrose_production(
        s.carbohydrates.starch_partition_coeff,
        s.photo.photosynthesis,
        s.carbohydrates.starch_degradation_rate,
        s.carbohydrates.uptake_cost,
        s.carbohydrates.transport_cost,
        s.carbohydrates.respiration_frequency,
        s.carbohydrates.sucrose,
        s.carbohydrates.sucrose_loading_frequency,
        s.carbohydrates.night_efficiency_starch);

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
}

void simulation_plant_rhs_adapter(
    real_t t,
    const void* state,
    void* dst,
    void* ctx)
{
    const real_t* x = (const real_t*)state;
    real_t* dxdt = (real_t*)dst;

    simulation_plant_rhs(t, x, dxdt, ctx);
}

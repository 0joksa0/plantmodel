#include "simulation/dynamics.h"

#include "model/model.h"

#include <math.h>

void simulation_pack_state(real_t* x, const Input* in)
{
    x[X_STARCH_PART] = in->core.starch_partition_coeff;
    x[X_STARCH] = in->core.starch;
    x[X_SUCROSE] = in->core.sucrose;
    x[X_N_AFF] = in->core.nitrogen_affinity;
    x[X_P_AFF] = in->core.phosphorus_affinity;
    x[X_N] = in->core.nitrogen;
    x[X_P] = in->core.phosphorus;
    x[X_SUCROSE_ROOT_ALLOC] = in->core.sucrose_root_allocation;
    x[X_LEAF_BIOMASS] = in->core.leaf_biomass;
    x[X_ROOT_BIOMASS] = in->core.root_biomass;
}

void simulation_unpack_state(Input* in, const real_t* x)
{
    in->core.starch_partition_coeff = x[X_STARCH_PART];
    in->core.starch = x[X_STARCH];
    in->core.sucrose = x[X_SUCROSE];
    in->core.nitrogen_affinity = x[X_N_AFF];
    in->core.phosphorus_affinity = x[X_P_AFF];
    in->core.nitrogen = x[X_N];
    in->core.phosphorus = x[X_P];
    in->core.sucrose_root_allocation = x[X_SUCROSE_ROOT_ALLOC];
    in->core.leaf_biomass = x[X_LEAF_BIOMASS];
    in->core.root_biomass = x[X_ROOT_BIOMASS];
    in->core.total_biomass = in->core.leaf_biomass + in->core.root_biomass;
}

void simulation_compute_algebraic(
    Input* s,
    real_t hour_of_day,
    real_t starch_night_start)
{
    update_light_conditions(s, hour_of_day);

    s->core.min_nitrogen = min_nitrogen(
        s->core.respiration_frequency, s->core.max_sucrose, s->core.max_starch,
        s->core.sucrose_loading_frequency, s->core.assimilation_cost_nitrogen,
        s->core.min_nitrogen_photosynthesis, s->core.nutrient_conversion_parameter, s->core.photoperiod);

    s->core.min_phosphorus = min_phosphorus(
        s->core.respiration_frequency, s->core.max_sucrose, s->core.max_starch,
        s->core.sucrose_loading_frequency, s->core.assimilation_cost_phosphorus,
        s->core.min_phosphorus_photosynthesis, s->core.nutrient_conversion_parameter, s->core.photoperiod);

    s->core.max_starch = max_starch(s->core.max_starch_degradation_rate, s->core.photoperiod);
    s->core.max_nitrogen = max_nitrogen(s->core.min_nitrogen);
    s->core.max_phosphorus = max_phosphorus(s->core.max_nitrogen, s->core.optimal_stoichiometric_ratio);

    s->core.nitrogen_saturation = nitrogen_saturation(s->core.nitrogen, s->core.min_nitrogen_photosynthesis);
    s->core.phosphorus_saturation = phosphorus_saturation(
        s->core.min_nitrogen_photosynthesis, s->core.optimal_stoichiometric_ratio, s->core.phosphorus);

    s->core.limitation_of_photosynthetic_rate = limitation_of_photosynthetic_rate(s->core.starch, s->core.feedback_on_photosynthesis, s->core.max_starch);

    s->core.photosynthesis = farquhar_photosynthesis(s);
    // s->core.photosynthesis = photosynthesis(s->core.light, s->core.limitation_of_photosynthetic_rate, s->core.max_photosynthetic_rate, s->core.nitrogen_saturation, s->core.phosphorus_saturation, s->core.leaf_biomass, s->core.min_leaf_biomass);

    s->core.night_efficiency_starch = night_efficiency_starch(
        s->core.sucrose, s->core.max_sucrose, s->core.lambda_g, s->core.light, s->core.starch_partition_coeff);

    s->core.starch_degradation_rate = starch_degradation(
        s->core.max_starch_degradation_rate,
        s->core.max_sucrose,
        starch_night_start,
        s->core.min_starch,
        s->core.sucrose,
        s->core.light,
        hour_of_day,
        s->core.photoperiod);

    s->core.uptake_cost = uptake_cost(
        s->core.nitrogen_uptake_sucrose_consumption, s->core.nitrogen_uptake,
        s->core.phosphorus_uptake_sucrose_consumption, s->core.phosphorus_uptake);

    s->core.transport_cost = transport_cost(
        s->core.sucrose_consumption_transport, s->core.respiration_frequency,
        s->core.sucrose, s->core.sucrose_loading_frequency, s->core.night_efficiency_starch);

    s->core.nitrogen_cost = nitrogen_cost(
        s->core.respiration_frequency, s->core.sucrose, s->core.starch,
        s->core.sucrose_loading_frequency, s->core.night_efficiency_starch,
        s->core.assimilation_cost_nitrogen,
        s->core.leaf_biomass, s->core.total_biomass,
        s->core.photosynthesis, s->core.max_photosynthetic_rate,
        s->core.min_nitrogen_photosynthesis);

    s->core.phosphorus_cost = phosphorus_cost(
        s->core.respiration_frequency, s->core.sucrose, s->core.starch,
        s->core.sucrose_loading_frequency, s->core.night_efficiency_starch,
        s->core.assimilation_cost_phosphorus,
        s->core.leaf_biomass, s->core.total_biomass,
        s->core.photosynthesis, s->core.max_photosynthetic_rate,
        s->core.min_phosphorus_photosynthesis);

    s->core.nitrogen_nutrient_uptake = nitrogen_nutrient_uptake(
        s->core.max_nitrogen_uptake, s->core.nitrogen_soil_content, s->core.Michaelis_Menten_constant_nitrogen);

    s->core.phosphorus_nutrient_uptake = phosphorus_nutrient_uptake(
        s->core.max_phosphorus_uptake, s->core.phosphorus_soil_content, s->core.Michaelis_Menten_constant_phosphorus);

    s->core.pot_nitrogen_uptake = pot_nitrogen_uptake(
        s->core.nitrogen_nutrient_uptake, s->core.nitrogen_affinity, s->core.root_biomass, s->core.total_biomass);

    s->core.pot_phosphorus_uptake = pot_phosphorus_uptake(
        s->core.phosphorus_nutrient_uptake, s->core.phosphorus_affinity, s->core.root_biomass, s->core.total_biomass);

    s->core.nitrogen_uptake = nitrogen_uptake(
        s->core.pot_nitrogen_uptake, s->core.sucrose, s->core.nitrogen_uptake_sucrose_consumption);

    s->core.phosphorus_uptake = phosphorus_uptake(
        s->core.pot_phosphorus_uptake, s->core.sucrose, s->core.phosphorus_uptake_sucrose_consumption);

    s->core.stoichiometric_signal = stoichiometric_signal(
        s->core.optimal_stoichiometric_ratio, s->core.nitrogen, s->core.phosphorus);
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
    const void* state,
    void* dst,
    void* ctx)
{
    const real_t* x = (const real_t*)state;
    real_t* dxdt = (real_t*)dst;

    simulation_plant_rhs(t, x, dxdt, ctx);
}

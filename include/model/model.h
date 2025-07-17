#ifndef MODEL_H
#define MODEL_H
#include "model/input.h"
#include <solver.h>

real_t photosynthesis(
    real_t light,
    real_t limitation_of_photosyntethic_rate,
    real_t max_photosyntetic_rate,
    real_t nitrogen_saturation,
    real_t phosphorus_saturation,
    real_t leaf_biomass,
    real_t min_leaf_biomass);
real_t farquhar_photosynthesis(Input* input);

real_t nitrogen_saturation(
    real_t nitrogen,
    real_t min_nitrogen_photosynthesis);

real_t phosphorus_saturation(
    real_t min_nitrogen_photosynthesis,
    real_t optimal_stechimetric_ratio,
    real_t phosphorus);

real_t limitation_of_photosyntethic_rate(
    real_t starch,
    real_t feedback_on_photosynthesis,
    real_t max_starch);

real_t max_starch(
    real_t max_starch_degradation_rate,
    real_t photoperiod);

real_t starch_degradation(
    real_t max_starch_degradation_rate,
    real_t max_sucrose,
    real_t starch_night_start,
    real_t min_starch,
    real_t sucrose,
    real_t light,
    real_t time,
    real_t photoperiod);

real_t starch_sucrose_partition(
    real_t light,
    real_t starch_partition_coeff,
    real_t min_sucrose,
    real_t max_sucrose,
    real_t sucrose,
    real_t lambda_sdr,
    real_t lambda_sdi,
    real_t lambda_sni);

real_t starch_sucrose_partition_f(
    real_t t,
    real_t starch_partition_coeff,
    void* params);

real_t nitrogen_content(
    real_t nitrogen_uptake,
    real_t respiration_frequency,
    real_t sucrose,
    real_t starch,
    real_t sucrose_loading_frequency,
    real_t night_efficiency_starch,
    real_t assimilaton_cost,
    real_t leaf_biomass,
    real_t total_biomass,
    real_t photosyntesis,
    real_t max_photosyntetic_rate,
    real_t nutrient_conversion_parameter,
    real_t min_nitrogen_photosynthesis);
real_t nitrogen_content_f(
    real_t t,
    real_t nitrogen,
    void* params);
real_t phosphorus_content(
    real_t phosphorus_uptake,
    real_t respiration_frequency,
    real_t sucrose,
    real_t starch,
    real_t sucrose_loading_frequency,
    real_t night_efficiency_starch,
    real_t assimilaton_cost,
    real_t leaf_biomass,
    real_t total_biomass,
    real_t photosyntesis,
    real_t max_photosyntetic_rate,
    real_t nutrient_conversion_parameter,
    real_t min_phosphorus_photosynthesis);
real_t phosphorus_content_f(
    real_t t,
    real_t phosphorus,
    void* params);

real_t assimilation_cost_nitrogen(
    real_t lambda_csn);

real_t assimilation_cost_phosphorus(
    real_t assimilation_cost_nitrogen,
    real_t optimal_stechiometric_ratio);

real_t pot_nitrogen_uptake(
    real_t nitrogen_nutrient_uptake,
    real_t nitrogen_affinity,
    real_t root_biomass,
    real_t total_biomass);

real_t nitrogen_uptake(
    real_t pot_nitrogen_uptake,
    real_t sucrose,
    real_t nitrogen_uptake_sucrose_consumption);

real_t nitrogen_nutrient_uptake(
    real_t max_nitrogen_uptake,
    real_t nitrogen_soil_content,
    real_t Michaelis_Menten_const_nitrogen);

real_t pot_phosphorus_uptake(
    real_t phosphorus_nutrient_uptake,
    real_t phosphorus_affinity,
    real_t root_biomass,
    real_t total_biomass);

real_t phosphorus_uptake(
    real_t pot_phosphorus_uptake,
    real_t sucrose,
    real_t phosphorus_uptake_sucrose_consumption);

real_t phosphorus_nutrient_uptake(
    real_t max_phosphorus_uptake,
    real_t phosphorus_soil_content,
    real_t Michaelis_Menten_const_phosphorus);
real_t nitrogen_affinity(
    real_t nitrogen_affinity,
    real_t nitrogen_uptake,
    real_t nitrogen_cost,
    real_t max_nitrogen,
    real_t nitrogen,
    real_t photosynthesis,
    real_t max_photosyntetic_rate,
    real_t lambda_k,
    real_t optimal_stochiometric_ratio,
    real_t phosphorus,
    real_t min_nitrogen,
    real_t nitrogen_uptake_sucrose_consumption,
    real_t starch_partition_coeff,
    real_t starch_degradation_rate);
real_t nitrogen_affinity_f(
    real_t t,
    real_t nitrogen_affinity_p,
    void* params);
real_t phosphorus_affinity(
    real_t phosphorus_affinity,
    real_t phosphorus_uptake,
    real_t phosphorus_cost,
    real_t max_phosphorus,
    real_t phosphorus,
    real_t photosynthesis,
    real_t max_photosyntetic_rate,
    real_t lambda_k,
    real_t optimal_stochiometric_ratio,
    real_t nitrogen,
    real_t min_phosphorus,
    real_t phosphorus_uptake_sucrose_consumption,
    real_t starch_partition_coeff,
    real_t starch_degradation_rate);
real_t phosphorus_affinity_f(
    real_t t,
    real_t phosphorus_affinity_p,
    void* params);
real_t nitrogen_cost(
    real_t respiration_frequency,
    real_t sucrose,
    real_t starch,
    real_t sucrose_loading_frequency,
    real_t night_efficiency_starch,
    real_t assimilation_cost_nitrogen,
    real_t leaf_biomass,
    real_t total_biomass,
    real_t photosynthesis,
    real_t max_photosyntetic_rate,
    real_t min_nitrogen_photosyntesis);

real_t phosphorus_cost(
    real_t respiration_frequency,
    real_t sucrose,
    real_t starch,
    real_t sucrose_loading_frequency,
    real_t night_efficiency_starch,
    real_t assimilation_cost_phosphorus,
    real_t leaf_biomass,
    real_t total_biomass,
    real_t photosynthesis,
    real_t max_photosyntetic_rate,
    real_t min_phosphorus_photosyntesis);

real_t min_nitrogen(
    real_t respiration_frequency,
    real_t max_sucrose,
    real_t max_starch,
    real_t sucrose_loading_frequency,
    real_t assimilation_cost_nitrogen,
    real_t min_nitrogen_photosynthesis,
    real_t nutrient_conversion_parameter,
    real_t photoperiod);

real_t min_phosphorus(
    real_t respiration_frequency,
    real_t max_sucrose,
    real_t max_starch,
    real_t sucrose_loading_frequency,
    real_t assimilation_cost_phosphorus,
    real_t min_phosphorus_photosynthesis,
    real_t nutrient_conversion_parameter,
    real_t photoperiod);

real_t max_nitrogen(real_t min_nitrogen);

real_t max_phosphorus(real_t max_nitrogen, real_t optimal_stechiometric_ratio);

real_t starch_production(
    real_t photosynthesis,
    real_t starch_partition_coeff,
    real_t starch_degradation_rate);

real_t starch_production_f(
    real_t t,
    real_t starch,
    void* params);

real_t sucrose_production(
    real_t sucrose_partition_coeff,
    real_t photosynthesis,
    real_t starch_degradation_rate,
    real_t uptake_cost,
    real_t transport_cost,
    real_t respiration_frequency,
    real_t sucrose,
    real_t sucrose_loading_frequency,
    real_t night_efficiency_starch);

real_t sucrose_production_f(
    real_t t,
    real_t sucrose,
    void* params);

real_t night_efficieny_starch(
    real_t sucrose,
    real_t max_sucrose,
    real_t lambda_g,
    real_t light,
    real_t starch_partition_coeff);

real_t uptake_cost(
    real_t nitrogen_uptake_sucrose_consumption,
    real_t nitrogen_uptake,
    real_t phosphorus_uptake_sucrose_consumption,
    real_t phosphorus_uptake);

real_t transport_cost(
    real_t sucrose_consumption_transport,
    real_t respiration_frequency,
    real_t sucrose,
    real_t sucrose_loading_frequency,
    real_t night_efficiency_starch);

real_t leaf_growth(
    real_t labmda_sb,
    real_t sucrose_root_allocation,
    real_t sucrose_loading_frequency,
    real_t night_efficiency_starch,
    real_t leaf_biomass,
    real_t leaf_death_rate,
    real_t leaf_competative_rate);
real_t leaf_growth_f(
    real_t t,
    real_t leaf_biomass,
    void* params);
real_t root_growth(
    real_t labmda_sb,
    real_t sucrose_root_allocation,
    real_t sucrose_loading_frequency,
    real_t night_efficiency_starch,
    real_t leaf_biomass,
    real_t root_biomass,
    real_t root_death_rate,
    real_t root_competative_rate);
real_t root_growth_f(
    real_t t,
    real_t root_biomass,
    void* params);
real_t stochiometric_signal(
    real_t optimal_stechiometric_ratio,
    real_t nitrogen,
    real_t phosphorus);

// resource_allocation -> tissue_priority
real_t sucrose_root_allocation(
    real_t sucrose_root_allocation,
    real_t nitrogen_affinity,
    real_t phospsorus_affinity,
    real_t stochiometric_signal,
    real_t sucrose,
    real_t min_sucrose,
    real_t nitrogen,
    real_t min_nitrogen,
    real_t phosphorus,
    real_t min_phosphorus);
real_t sucrose_root_allocation_f(
    real_t t,
    real_t sucrose_root_allocation_p,
    void* params);
// Farquhar C3 Photosynthesis Model

real_t rubisco_limited_photosynthesis(
    real_t intercellular_CO2, // Ci
    real_t CO2_compensation_point, // Gamma_star
    real_t max_rubisco_carboxylation_rate, // Vcmax
    real_t michaelis_constant_CO2, // Kc
    real_t michaelis_constant_O2, // Ko
    real_t oxygen_concentration // O
);

// Light-limited (RuBP regeneration) assimilation rate
real_t light_limited_photosynthesis(
    real_t intercellular_CO2, // Ci
    real_t CO2_compensation_point, // Gamma_star
    real_t electron_transport_rate // J
);

// Net photosynthesis rate
real_t net_photosynthesis(
    real_t rubisco_limited_rate, // Ac
    real_t light_limited_rate, // Aj
    real_t respiration_rate // Rd
);

// Convert An (umol CO2 m^-2 s^-1) to ph (umol C6 gFW^-1 h^-1)
real_t convert_to_photosynthesis_per_gram_per_hour(
    real_t net_photosynthesis_rate, // An [μmol CO2 / m² / s]
    real_t leaf_biomass, // gFW
    real_t specific_leaf_area // SLA [m² / gFW]
);

// Intercellular CO2 based on diffusion from ambient and net photosynthesis
real_t intercellular_CO2(
    real_t ambient_CO2_concentration, // Ca
    real_t net_photosynthesis_rate, // An
    real_t stomatal_conductance // gs
);

// Maksimov
real_t moistureAvailability(
    real_t evapotranspiration,
    real_t slope_runoff,
    real_t vertial_filtration_down,
    real_t vertial_filtration_up,
    real_t initial_soil_moisture_reserve,
    real_t precipitation);

real_t lightAvailability(
    real_t available_light_coef,
    real_t leaf_radiation,
    real_t soil_radiation);

real_t heatAvailability(
    real_t heat_capacity_air_layer,
    real_t temp_change_air_layer,
    real_t heat_capacity_plant_mass,
    real_t temp_change_plant_layer,
    real_t heat_capacity_root_zone_soil,
    real_t temp_change_root_zone,
    real_t thermal_conductivity_soil,
    real_t temp_gradient_deep_soil,
    real_t shortwave_radiation,
    real_t longwave_radiation);

real_t foodSecurity(
    real_t sum_total_nutrients_on_slope,
    real_t sum_mastered_by_system,
    real_t sum_initial_nutrients,
    real_t sum_applied_nutrients,
    real_t sum_hardened_nutrients,
    real_t sum_lost_nutrients,
    int is_accumulation_zone);

real_t compute_gas_supply_coefficient(
    real_t CO2_density,
    real_t turbulent_diffusion_coefficient,
    real_t dCO2_soil_dx,
    real_t dCO2_air_dx,
    real_t co2_production_plant_layer);

real_t growthPotential(
    real_t power_spent,
    real_t mass);

real_t iterate_ci(Input* input, int max_iter, real_t epsilon);
#endif

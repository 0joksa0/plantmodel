#ifndef MODEL_H
#define MODEL_H

typedef struct {
    float light;
    float limitation_of_photosyntetic_rate;
    float max_photosyntetic_rate;
    float nitrogen_saturation;
    float phosphorus_saturation;
    float leaf_biomass;
    float min_leaf_biomass;
    float feedback_on_photosynthesis;
    float max_starch;
    float max_starch_degradation_rate;
    float photoperiod;
    float optimal_stechiometric_ratio;
    float photosynthesis;
    float starch_partition_coeff;
    float starch_degradation_rate;
    float uptake_cost;
    float transport_cost;
    float respiration_frequency;
    float sucrose_loading_frequency;
    float night_efficiency_starch;
    float max_sucrose;
    float min_sucrose;
    float starch_night_start;
    float min_starch;
    float lambda_sdr;
    float lambda_sdi;
    float lambda_sni;
    float lambda_g;
    float lambda_sb;
    float sucrose;
    float starch;
    float nitrogen;
    float phosphorus;
    float nitrogen_uptake;
    float phosphorus_uptake;
    float assimilation_cost_nitrogen;
    float lambda_csn;
    float min_nitrogen_photosynthesis;
    float min_phosphorus_photosynthesis;
    float nutrient_conversion_parameter;
    float total_biomass;
    float assimilation_cost_phorphorus;
    float pot_nitrogen_uptake;
    float pot_phosphorus_uptake;
    float nitrogen_nutrient_uptake;
    float nitrogen_affinity;
    float phosphorus_nutrient_uptake;
    float phosphorus_affinity;
    float root_biomass;
    float nitrogen_uptake_sucrose_consumption;
    float phosphorus_uptake_sucrose_consumption;
    float max_nitrogen_uptake;
    float max_phosphorus_uptake;
    float Michaelis_Menten_constant_nitrogen;
    float Michaelis_Menten_constant_phosphorus;
    float nitrogen_soil_content;
    float phosphorus_soil_content;
    float lambda_k;
    float max_nitrogen;
    float min_nitrogen;
    float max_phosphorus;
    float min_phosphorus;
    float nitrogen_cost;
    float phosphorus_cost;
    float sucrose_consumption_transport;
    float sucrose_root_allocation;
    float stochiometric_signal;
    float leaf_deathrate;
    float root_deathrate;
    float leaf_competitive_rate;
    float root_competitive_rate;

} Input;

float photosynthesis(
    float light,
    float limitation_of_photosyntethic_rate,
    float max_photosyntetic_rate,
    float nitrogen_saturation,
    float phosphorus_saturation,
    float leaf_biomass,
    float min_leaf_biomass);

float nitrogen_saturation(
    float nitrogen,
    float min_nitrogen_photosynthesis);

float phosphorus_saturation(
    float min_nitrogen_photosynthesis,
    float optimal_stechimetric_ratio,
    float phosphorus);

float limitation_of_photosyntethic_rate(
    float starch,
    float feedback_on_photosynthesis,
    float max_starch);

float max_starch(
    float max_starch_degradation_rate,
    float photoperiod);

float starch_degradation(
    float max_starch_degradation_rate,
    float max_sucrose,
    float starch_night_start,
    float min_starch,
    float sucrose,
    float light,
    float time,
    float photoperiod);

float starch_sucrose_partition(
    float light,
    float starch_partition_coeff,
    float min_sucrose,
    float max_sucrose,
    float sucrose,
    float lambda_sdr,
    float lambda_sdi,
    float lambda_sni);

float nitrogen_content(
    float nitrogen_uptake,
    float respiration_frequency,
    float sucrose,
    float starch,
    float sucrose_loading_frequency,
    float night_efficiency_starch,
    float assimilaton_cost,
    float leaf_biomass,
    float total_biomass,
    float photosyntesis,
    float max_photosyntetic_rate,
    float nutrient_conversion_parameter,
    float min_nitrogen_photosynthesis);

float phosphorus_content(
    float phosphorus_uptake,
    float respiration_frequency,
    float sucrose,
    float starch,
    float sucrose_loading_frequency,
    float night_efficiency_starch,
    float assimilaton_cost,
    float leaf_biomass,
    float total_biomass,
    float photosyntesis,
    float max_photosyntetic_rate,
    float nutrient_conversion_parameter,
    float min_phosphorus_photosynthesis);

float assimilation_cost_nitrogen(
    float lambda_csn);

float assimilation_cost_phosphorus(
    float assimilation_cost_nitrogen,
    float optimal_stechiometric_ratio);

float pot_nitrogen_uptake(
    float nitrogen_nutrient_uptake,
    float nitrogen_affinity,
    float root_biomass,
    float total_biomass);

float nitrogen_uptake(
    float pot_nitrogen_uptake,
    float sucrose,
    float nitrogen_uptake_sucrose_consumption);

float nitrogen_nutrient_uptake(
    float max_nitrogen_uptake,
    float nitrogen_soil_content,
    float Michaelis_Menten_const_nitrogen);

float pot_phosphorus_uptake(
    float phosphorus_nutrient_uptake,
    float phosphorus_affinity,
    float root_biomass,
    float total_biomass);

float phosphorus_uptake(
    float pot_phosphorus_uptake,
    float sucrose,
    float phosphorus_uptake_sucrose_consumption);

float phosphorus_nutrient_uptake(
    float max_phosphorus_uptake,
    float phosphorus_soil_content,
    float Michaelis_Menten_const_phosphorus);
float nitrogen_affinity(
    float nitrogen_affinity,
    float nitrogen_uptake,
    float nitrogen_cost,
    float max_nitrogen,
    float nitrogen,
    float photosynthesis,
    float max_photosyntetic_rate,
    float lambda_k,
    float optimal_stochiometric_ratio,
    float phosphorus,
    float min_nitrogen,
    float nitrogen_uptake_sucrose_consumption,
    float starch_partition_coeff,
    float starch_degradation_rate);

float phosphorus_affinity(
    float phosphorus_affinity,
    float phosphorus_uptake,
    float phosphorus_cost,
    float max_phosphorus,
    float phosphorus,
    float photosynthesis,
    float max_photosyntetic_rate,
    float lambda_k,
    float optimal_stochiometric_ratio,
    float nitrogen,
    float min_phosphorus,
    float phosphorus_uptake_sucrose_consumption,
    float starch_partition_coeff,
    float starch_degradation_rate);

float nitrogen_cost(
    float respiration_frequency,
    float sucrose,
    float starch,
    float sucrose_loading_frequency,
    float night_efficiency_starch,
    float assimilation_cost_nitrogen,
    float leaf_biomass,
    float total_biomass,
    float photosynthesis,
    float max_photosyntetic_rate,
    float min_nitrogen_photosyntesis);

float phosphorus_cost(
    float respiration_frequency,
    float sucrose,
    float starch,
    float sucrose_loading_frequency,
    float night_efficiency_starch,
    float assimilation_cost_phosphorus,
    float leaf_biomass,
    float total_biomass,
    float photosynthesis,
    float max_photosyntetic_rate,
    float min_phosphorus_photosyntesis);

float min_nitrogen(
    float respiration_frequency,
    float max_sucrose,
    float max_starch,
    float sucrose_loading_frequency,
    float assimilation_cost_nitrogen,
    float min_nitrogen_photosynthesis,
    float nutrient_conversion_parameter,
    float photoperiod);

float min_phosphorus(
    float respiration_frequency,
    float max_sucrose,
    float max_starch,
    float sucrose_loading_frequency,
    float assimilation_cost_phosphorus,
    float min_phosphorus_photosynthesis,
    float nutrient_conversion_parameter,
    float photoperiod);

float max_nitrogen(float min_nitrogen);

float max_phosphorus(float max_nitrogen, float optimal_stechiometric_ratio);

float starch_production(
    float photosynthesis,
    float starch_partition_coeff,
    float starch_degradation_rate);

float sucrose_production(
    float sucrose_partition_coeff,
    float photosynthesis,
    float starch_degradation_rate,
    float uptake_cost,
    float transport_cost,
    float respiration_frequency,
    float sucrose,
    float sucrose_loading_frequency,
    float night_efficiency_starch);

float night_efficieny_starch(
    float sucrose,
    float max_sucrose,
    float lambda_g,
    float light,
    float starch_partition_coeff);

float uptake_cost(
    float nitrogen_uptake_sucrose_consumption,
    float nitrogen_uptake,
    float phosphorus_uptake_sucrose_consumption,
    float phosphorus_uptake);

float transport_cost(
    float sucrose_consumption_transport,
    float respiration_frequency,
    float sucrose,
    float sucrose_loading_frequency,
    float night_efficiency_starch);

float leaf_growth(
    float labmda_sb,
    float sucrose_root_allocation,
    float sucrose_loading_frequency,
    float night_efficiency_starch,
    float leaf_biomass,
    float leaf_death_rate,
    float leaf_competative_rate);

float root_growth(
    float labmda_sb,
    float sucrose_root_allocation,
    float sucrose_loading_frequency,
    float night_efficiency_starch,
    float leaf_biomass,
    float root_biomass,
    float root_death_rate,
    float root_competative_rate);

float stochiometric_signal(
    float optimal_stechiometric_ratio,
    float nitrogen,
    float phosphorus);

// resource_allocation -> tissue_priority
float sucrose_root_allocation(
    float sucrose_root_allocation,
    float nitrogen_affinity,
    float phospsorus_affinity,
    float stochiometric_signal,
    float sucrose,
    float min_sucrose,
    float nitrogen,
    float min_nitrogen,
    float phosphorus,
    float min_phosphorus);

// Maksimov
float moistureAvailability(
    float evapotranspiration,
    float slope_runoff,
    float vertial_filtration_down,
    float vertial_filtration_up,
    float initial_soil_moisture_reserve,
    float precipitation);

float lightAvailability(
    float available_light_coef,
    float leaf_radiation,
    float soil_radiation);

float heatAvailability(
    float heat_capacity_air_layer,
    float temp_change_air_layer,
    float heat_capacity_plant_mass,
    float temp_change_plant_layer,
    float heat_capacity_root_zone_soil,
    float temp_change_root_zone,
    float thermal_conductivity_soil,
    float temp_gradient_deep_soil,
    float shortwave_radiation,
    float longwave_radiation);

float foodSecurity(
    float sum_total_nutrients_on_slope,
    float sum_mastered_by_system,
    float sum_initial_nutrients,
    float sum_applied_nutrients,
    float sum_hardened_nutrients,
    float sum_lost_nutrients,
    int is_accumulation_zone);

float compute_gas_supply_coefficient(
    float CO2_density,
    float turbulent_diffusion_coefficient,
    float dCO2_soil_dx,
    float dCO2_air_dx,
    float co2_production_plant_layer);

float growthPotential(
    float power_spent,
    float mass);

#endif

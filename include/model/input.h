#ifndef INPUT_H
#define INPUT_H


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


Input generate_input();

#endif

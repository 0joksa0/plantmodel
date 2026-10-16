#ifndef INPUT_H
#define INPUT_H
#include <solver.h>

#define TEMPERATURE_UNSET REAL(-900.0)

typedef struct {
    real_t light;
    real_t light_PAR;
    real_t max_PAR;
    real_t light_profile;  /* 0 = half-sine (default), 1 = square */
    real_t limitation_of_photosynthetic_rate;
    real_t max_photosynthetic_rate;
    real_t nitrogen_saturation;
    real_t phosphorus_saturation;
    real_t min_leaf_biomass;
    real_t extinction_coefficient;  /* canopy light extinction [-] */
    real_t feedback_on_photosynthesis;
    real_t photoperiod;
    real_t photosynthesis;
    int    photo_mode;
} InputPhoto;

typedef struct {
    real_t max_starch;
    real_t max_starch_degradation_rate;
    real_t starch_partition_coeff;
    real_t starch_degradation_rate;
    real_t uptake_cost;
    real_t transport_cost;
    real_t respiration_frequency;
    real_t sucrose_loading_frequency;
    real_t night_efficiency_starch;
    real_t sink_activation_half;   /* half-activation sucrose [umolC6/gFW] */
    real_t sink_activation_hill;   /* Hill coefficient [-] */
    real_t max_sucrose;
    real_t min_sucrose;
    real_t starch_night_start;
    real_t min_starch;
    real_t lambda_sdr;
    real_t lambda_sdi;
    real_t lambda_sni;
    real_t lambda_g;
    real_t sucrose;
    real_t starch;
    real_t sucrose_consumption_transport;

    /* Vacuolar sucrose compartment.  Whole-tissue sucrose measurements
     * (Sulpice et al. 2014; Annunziata et al. 2017) pool cytosol and vacuole,
     * while the model tracked the cytosolic pool alone.  The partition ratio is
     * derived from Seydel et al. non-aqueous fractionation; see
     * config/parameter_provenance.csv for the tiering of both rate constants. */
    real_t sucrose_vacuole;  /* vacuolar sucrose [umol/gFW]; state variable */
    real_t k_vacuole_in;     /* cytosol -> vacuole rate constant [h-1] */
    real_t k_vacuole_out;    /* vacuole -> cytosol rate constant [h-1] */
} InputCarbohydrates;

typedef struct {
    real_t optimal_stoichiometric_ratio;
    real_t nitrogen;
    real_t phosphorus;
    real_t nitrogen_uptake;
    real_t phosphorus_uptake;
    real_t assimilation_cost_nitrogen;
    real_t lambda_csn;
    real_t min_nitrogen_photosynthesis;
    real_t min_phosphorus_photosynthesis;
    real_t nutrient_conversion_parameter;
    real_t assimilation_cost_phosphorus;
    real_t pot_nitrogen_uptake;
    real_t pot_phosphorus_uptake;
    real_t nitrogen_nutrient_uptake;
    real_t nitrogen_affinity;
    real_t phosphorus_nutrient_uptake;
    real_t phosphorus_affinity;
    real_t nitrogen_uptake_sucrose_consumption;
    real_t phosphorus_uptake_sucrose_consumption;
    real_t max_nitrogen_uptake;
    real_t max_phosphorus_uptake;
    real_t Michaelis_Menten_constant_nitrogen;
    real_t Michaelis_Menten_constant_phosphorus;
    real_t nitrogen_soil_content;
    real_t phosphorus_soil_content;
    real_t lambda_k;
    real_t max_nitrogen;
    real_t min_nitrogen;
    real_t max_phosphorus;
    real_t min_phosphorus;
    real_t nitrogen_cost;
    real_t phosphorus_cost;
    real_t stoichiometric_signal;
} InputNutrients;

typedef struct {
    real_t leaf_biomass;
    real_t root_biomass;
    real_t total_biomass;
    real_t lambda_sb;
    real_t sucrose_root_allocation;
    real_t leaf_deathrate;
    real_t root_deathrate;
    real_t leaf_competitive_rate;
    real_t root_competitive_rate;

    /* Construction cost of new biomass, in umol hexose per gram fresh
     * weight built.  Y_G / (M_hexose * 1e-6 * FW:DW) with Y_G = 1.4 g
     * glucose gDW-1 (Penning de Vries et al., growth respiration) and
     * FW:DW = 12; see growth_cost() in src/model/carbohydrates.c.
     * A value of 0 reproduces the previous uncharged behaviour. */
    real_t growth_stoichiometry;
} InputGrowth;

typedef struct {
    real_t intercellular_CO2;
    real_t ambient_CO2_concentration;
    real_t net_photosynthesis_rate;
    real_t stomatal_conductance;
    real_t convert_to_photosynthesis_per_gram_per_hour;
    real_t specific_leaf_area;
    real_t net_photosynthesis;
    real_t rubisco_limited_rate;
    real_t light_limited_rate;
    real_t respiration_rate;
    real_t light_limited_photosynthesis;
    real_t CO2_compensation_point;
    real_t electron_transport_rate;
    real_t rubisco_limited_photosynthesis;
    real_t max_rubisco_carboxylation_rate;
    real_t michaelis_constant_CO2;
    real_t michaelis_constant_O2;
    real_t oxygen_concentration;
    real_t leaf_temperature;
    real_t leaf_temperature_day;    /* thermoperiod day set point [C] */
    real_t leaf_temperature_night;  /* thermoperiod night set point [C] */
    real_t relative_humidity;
    real_t leaf_water_potential;
} InputGasExchange;

typedef struct {
    InputPhoto photo;
    InputCarbohydrates carbohydrates;
    InputNutrients nutrients;
    InputGrowth growth;
    InputGasExchange gas_exchange;
} Input;

Input generate_input();

#endif

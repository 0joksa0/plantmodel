#include "model/model.h"

static real_t epsilon = REAL_EPSILON;

real_t nitrogen_saturation(
    real_t nitrogen,
    real_t min_nitrogen_photosynthesis)
{
    if (nitrogen < min_nitrogen_photosynthesis) {
        return 0;
    }
    return ((REAL(2.0) * nitrogen) / (nitrogen + min_nitrogen_photosynthesis)) - REAL(1.0);
}

real_t phosphorus_saturation(
    real_t min_nitrogen_photosynthesis,
    real_t optimal_stoichiometric_ratio,
    real_t phosphorus)
{
    real_t min_phosphorus_photosynthesis = min_nitrogen_photosynthesis / optimal_stoichiometric_ratio;
    if (phosphorus < min_phosphorus_photosynthesis) {
        return 0;
    }
    return ((REAL(2.0) * phosphorus) / (phosphorus + min_phosphorus_photosynthesis)) - REAL(1.0);
}

real_t nitrogen_content(
    real_t nitrogen_uptake,
    real_t respiration_frequency,
    real_t sucrose,
    real_t starch,
    real_t sucrose_loading_frequency,
    real_t night_efficiency_starch,
    real_t assimilation_cost,
    real_t leaf_biomass,
    real_t total_biomass,
    real_t photosynthesis,
    real_t max_photosynthetic_rate,
    real_t nutrient_conversion_parameter,
    real_t min_nitrogen_photosynthesis)
{

    return nitrogen_uptake - ((respiration_frequency * (sucrose + starch) + (night_efficiency_starch * sucrose_loading_frequency)) * assimilation_cost * (leaf_biomass / (total_biomass + epsilon))) - ((photosynthesis / max_photosynthetic_rate) * min_nitrogen_photosynthesis * nutrient_conversion_parameter);
}

real_t nitrogen_content_f(
    real_t t,
    real_t nitrogen,
    void* params)
{
    Input* input = (Input*)params;
    (void)t;
    (void)nitrogen;
    return nitrogen_content(
        input->nutrients.nitrogen_uptake,
        input->carbohydrates.respiration_frequency,
        input->carbohydrates.sucrose,
        input->carbohydrates.starch,
        input->carbohydrates.sucrose_loading_frequency,
        input->carbohydrates.night_efficiency_starch,
        input->nutrients.assimilation_cost_nitrogen,
        input->growth.leaf_biomass,
        input->growth.total_biomass,
        input->photo.photosynthesis,
        input->photo.max_photosynthetic_rate,
        input->nutrients.nutrient_conversion_parameter,
        input->nutrients.min_nitrogen_photosynthesis);
}

real_t phosphorus_content(
    real_t phosphorus_uptake,
    real_t respiration_frequency,
    real_t sucrose,
    real_t starch,
    real_t sucrose_loading_frequency,
    real_t night_efficiency_starch,
    real_t assimilation_cost,
    real_t leaf_biomass,
    real_t total_biomass,
    real_t photosynthesis,
    real_t max_photosynthetic_rate,
    real_t nutrient_conversion_parameter,
    real_t min_phosphorus_photosynthesis)
{

    return phosphorus_uptake - ((respiration_frequency * (sucrose + starch) + (night_efficiency_starch * sucrose_loading_frequency)) * assimilation_cost * (leaf_biomass / (total_biomass + epsilon))) - ((photosynthesis / max_photosynthetic_rate) * min_phosphorus_photosynthesis * nutrient_conversion_parameter);
}

real_t phosphorus_content_f(
    real_t t,
    real_t phosphorus,
    void* params)
{
    Input* input = (Input*)params;
    (void)t;
    (void)phosphorus;

    return phosphorus_content(
        input->nutrients.phosphorus_uptake,
        input->carbohydrates.respiration_frequency,
        input->carbohydrates.sucrose,
        input->carbohydrates.starch,
        input->carbohydrates.sucrose_loading_frequency,
        input->carbohydrates.night_efficiency_starch,
        input->nutrients.assimilation_cost_phosphorus,
        input->growth.leaf_biomass,
        input->growth.total_biomass,
        input->photo.photosynthesis,
        input->photo.max_photosynthetic_rate,
        input->nutrients.nutrient_conversion_parameter,
        input->nutrients.min_phosphorus_photosynthesis);
}

real_t assimilation_cost_nitrogen(
    real_t lambda_csn)
{
    return lambda_csn * REAL(1.724);
}

real_t assimilation_cost_phosphorus(
    real_t assimilation_cost_nitrogen,
    real_t optimal_stoichiometric_ratio)
{

    return assimilation_cost_nitrogen / optimal_stoichiometric_ratio;
}

real_t pot_nitrogen_uptake(
    real_t nitrogen_nutrient_uptake,
    real_t nitrogen_affinity,
    real_t root_biomass,
    real_t total_biomass)
{
    return nitrogen_nutrient_uptake * nitrogen_affinity * (root_biomass / (total_biomass + epsilon));
}

real_t nitrogen_uptake(
    real_t pot_nitrogen_uptake,
    real_t sucrose,
    real_t nitrogen_uptake_sucrose_consumption)
{

    return pot_nitrogen_uptake * (sucrose / (sucrose + (nitrogen_uptake_sucrose_consumption * pot_nitrogen_uptake * REAL(1.0)) + epsilon));
}

real_t nitrogen_nutrient_uptake(
    real_t max_nitrogen_uptake,
    real_t nitrogen_soil_content,
    real_t Michaelis_Menten_const_nitrogen)
{
    return max_nitrogen_uptake * (nitrogen_soil_content / (nitrogen_soil_content + Michaelis_Menten_const_nitrogen));
}

real_t pot_phosphorus_uptake(
    real_t phosphorus_nutrient_uptake,
    real_t phosphorus_affinity,
    real_t root_biomass,
    real_t total_biomass)
{
    return phosphorus_nutrient_uptake * phosphorus_affinity * (root_biomass / (total_biomass + epsilon));
}

real_t phosphorus_uptake(
    real_t pot_phosphorus_uptake,
    real_t sucrose,
    real_t phosphorus_uptake_sucrose_consumption)
{

    return pot_phosphorus_uptake * sucrose / (sucrose + (phosphorus_uptake_sucrose_consumption * pot_phosphorus_uptake * REAL(1.0)) + epsilon);
}

real_t phosphorus_nutrient_uptake(
    real_t max_phosphorus_uptake,
    real_t phosphorus_soil_content,
    real_t Michaelis_Menten_const_phosphorus)
{
    return max_phosphorus_uptake * (phosphorus_soil_content / (phosphorus_soil_content + Michaelis_Menten_const_phosphorus));
}

real_t nitrogen_affinity(
    real_t nitrogen_affinity,
    real_t nitrogen_uptake,
    real_t nitrogen_cost,
    real_t max_nitrogen,
    real_t nitrogen,
    real_t photosynthesis,
    real_t max_photosynthetic_rate,
    real_t lambda_k,
    real_t optimal_stoichiometric_ratio,
    real_t phosphorus,
    real_t min_nitrogen,
    real_t nitrogen_uptake_sucrose_consumption,
    real_t starch_partition_coeff,
    real_t starch_degradation_rate)
{
    real_t n_nmax = (max_nitrogen / (nitrogen + max_nitrogen));
    real_t un_real = (REAL(1.0) - (nitrogen_uptake / (nitrogen_uptake + nitrogen_cost + epsilon)));
    real_t n_op = nitrogen_affinity * lambda_k * ((nitrogen - (optimal_stoichiometric_ratio * phosphorus)) / (nitrogen + (optimal_stoichiometric_ratio * phosphorus) + epsilon));

    real_t term1 = (REAL(1.0) - nitrogen_affinity) * ((un_real * n_nmax) + (photosynthesis / max_photosynthetic_rate) - n_op);

    real_t n_nmin = (nitrogen / (nitrogen + min_nitrogen));
    real_t term2 = nitrogen_affinity * (n_nmin + ((nitrogen_uptake_sucrose_consumption * nitrogen_uptake) / ((nitrogen_uptake_sucrose_consumption * nitrogen_uptake) + (photosynthesis * (REAL(1.0) - starch_partition_coeff)) + starch_degradation_rate + epsilon)));
    return term1 - term2;
}

real_t nitrogen_affinity_f(
    real_t t,
    real_t nitrogen_affinity_p,
    void* params)
{
    Input* in = (Input*)params;
    (void)t;

    return nitrogen_affinity(nitrogen_affinity_p,
        in->nutrients.nitrogen_uptake,
        in->nutrients.nitrogen_cost,
        in->nutrients.max_nitrogen,
        in->nutrients.nitrogen,
        in->photo.photosynthesis,
        in->photo.max_photosynthetic_rate,
        in->nutrients.lambda_k,
        in->nutrients.optimal_stoichiometric_ratio,
        in->nutrients.phosphorus,
        in->nutrients.min_nitrogen,
        in->nutrients.nitrogen_uptake_sucrose_consumption,
        in->carbohydrates.starch_partition_coeff,
        in->carbohydrates.starch_degradation_rate);
}

real_t phosphorus_affinity(
    real_t phosphorus_affinity,
    real_t phosphorus_uptake,
    real_t phosphorus_cost,
    real_t max_phosphorus,
    real_t phosphorus,
    real_t photosynthesis,
    real_t max_photosynthetic_rate,
    real_t lambda_k,
    real_t optimal_stoichiometric_ratio,
    real_t nitrogen,
    real_t min_phosphorus,
    real_t phosphorus_uptake_sucrose_consumption,
    real_t starch_partition_coeff,
    real_t starch_degradation_rate)
{
    real_t n_nmax = (max_phosphorus / (phosphorus + max_phosphorus));
    real_t un_real = (REAL(1.0) - (phosphorus_uptake / (phosphorus_uptake + phosphorus_cost + epsilon)));
    real_t n_op = phosphorus_affinity * lambda_k * ((nitrogen - (optimal_stoichiometric_ratio * phosphorus)) / (nitrogen + (optimal_stoichiometric_ratio * phosphorus) + epsilon));

    real_t term1 = (REAL(1.0) - phosphorus_affinity) * ((un_real * n_nmax) + (photosynthesis / max_photosynthetic_rate) + n_op);

    real_t n_nmin = (phosphorus / (phosphorus + min_phosphorus));
    real_t term2 = phosphorus_affinity * (n_nmin + ((phosphorus_uptake_sucrose_consumption * phosphorus_uptake) / ((phosphorus_uptake_sucrose_consumption * phosphorus_uptake) + (photosynthesis * (REAL(1.0) - starch_partition_coeff)) + starch_degradation_rate + epsilon)));
    return term1 - term2;
}

real_t phosphorus_affinity_f(
    real_t t,
    real_t phosphorus_affinity_p,
    void* params)
{
    Input* input = (Input*)params;
    (void)t;

    return phosphorus_affinity(phosphorus_affinity_p, input->nutrients.phosphorus_uptake, input->nutrients.phosphorus_cost, input->nutrients.max_phosphorus, input->nutrients.phosphorus, input->photo.photosynthesis, input->photo.max_photosynthetic_rate, input->nutrients.lambda_k, input->nutrients.optimal_stoichiometric_ratio, input->nutrients.nitrogen, input->nutrients.min_phosphorus, input->nutrients.phosphorus_uptake_sucrose_consumption, input->carbohydrates.starch_partition_coeff, input->carbohydrates.starch_degradation_rate);
}

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
    real_t max_photosynthetic_rate,
    real_t min_nitrogen_photosynthesis)
{
    real_t first_term = (respiration_frequency * (sucrose + starch) + (night_efficiency_starch * sucrose_loading_frequency))
        * assimilation_cost_nitrogen * (leaf_biomass / (total_biomass + epsilon));

    real_t second_term = (photosynthesis / max_photosynthetic_rate) * min_nitrogen_photosynthesis;

    return first_term + second_term;
}

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
    real_t max_photosynthetic_rate,
    real_t min_phosphorus_photosynthesis)
{
    real_t first_term = (respiration_frequency * (sucrose + starch) + (night_efficiency_starch * sucrose_loading_frequency))
        * assimilation_cost_phosphorus * (leaf_biomass / (total_biomass + epsilon));

    real_t second_term = (photosynthesis / max_photosynthetic_rate) * min_phosphorus_photosynthesis;

    return first_term + second_term;
}

real_t min_nitrogen(
    real_t respiration_frequency,
    real_t max_sucrose,
    real_t max_starch,
    real_t sucrose_loading_frequency,
    real_t assimilation_cost_nitrogen,
    real_t min_nitrogen_photosynthesis,
    real_t nutrient_conversion_parameter,
    real_t photoperiod)
{

    real_t n_min = ((respiration_frequency * (max_sucrose + max_starch)) + (sucrose_loading_frequency * max_sucrose)) * (((REAL(1.0) + REAL(0.035)) * assimilation_cost_nitrogen) + (min_nitrogen_photosynthesis * nutrient_conversion_parameter * photoperiod));

    return n_min;
}

real_t min_phosphorus(
    real_t respiration_frequency,
    real_t max_sucrose,
    real_t max_starch,
    real_t sucrose_loading_frequency,
    real_t assimilation_cost_phosphorus,
    real_t min_phosphorus_photosynthesis,
    real_t nutrient_conversion_parameter,
    real_t photoperiod)
{

    real_t n_min = ((respiration_frequency * (max_sucrose + max_starch)) + (sucrose_loading_frequency * max_sucrose)) * (((REAL(1.0) + REAL(0.035)) * assimilation_cost_phosphorus) + (min_phosphorus_photosynthesis * nutrient_conversion_parameter * photoperiod));

    return n_min;
}

real_t max_nitrogen(real_t min_nitrogen)
{
    return min_nitrogen * REAL(4.0);
}

real_t max_phosphorus(real_t max_nitrogen, real_t optimal_stoichiometric_ratio)
{
    return max_nitrogen / optimal_stoichiometric_ratio;
}

real_t stoichiometric_signal(
    real_t optimal_stoichiometric_ratio,
    real_t nitrogen,
    real_t phosphorus)
{
    return optimal_stoichiometric_ratio / (optimal_stoichiometric_ratio + (nitrogen / (phosphorus + epsilon)));
}

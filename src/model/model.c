#include "model/model.h"
#include "model/input.h"
#include <math.h>
#include <solver.h>
#include <stdio.h>

real_t epsylon = REAL_EPSILON;

real_t photosynthesis(
    real_t light,
    real_t limitation_of_photosynthetic_rate,
    real_t max_photosyntetic_rate,
    real_t nitrogen_saturation,
    real_t phosphorus_saturation,
    real_t leaf_biomass,
    real_t min_leaf_biomass)
{
    if (leaf_biomass < min_leaf_biomass) {
        return 0;
    }

    // printf("lim_p_rate: %f, nit_sat: %f, phos_sat: %f\n", limitation_of_photosynthetic_rate, nitrogen_saturation, phosphorus_saturation);
    real_t photo = light * limitation_of_photosynthetic_rate * max_photosyntetic_rate * RMIN(nitrogen_saturation, phosphorus_saturation);
    return photo;
}

real_t farquhar_photosynthesis(Input* input)
{
    if (input->light <= REAL(0.0) || input->leaf_biomass < input->min_leaf_biomass) {
        return REAL(0.0);
    }

    real_t VPD = calculate_vapor_pressure_deficit(input->leaf_temperature, input->relative_humidity);
    input->stomatal_conductance = calculate_stomatal_conductance_jarvis(
        input->light_PAR,
        VPD,
        input->leaf_temperature,
        input->ambient_CO2_concentration,
        input->leaf_water_potential);
    // printf("stomal cond: %f ", input->stomatal_conductance);
    real_t An = iterate_ci(input, 50, REAL(0.0001));

    real_t limitation = input->limitation_of_photosyntetic_rate;
    real_t np_saturation = RMIN(input->nitrogen_saturation, input->phosphorus_saturation);
    // printf("lim:%f psat:%f nsat:%f ", input->limitation_of_photosyntetic_rate, input->phosphorus_saturation, input->nitrogen_saturation);

    input->convert_to_photosynthesis_per_gram_per_hour = convert_to_photosynthesis_per_gram_per_hour(
        An,
        input->leaf_biomass,
        input->specific_leaf_area);
    real_t adjusted_PH = input->convert_to_photosynthesis_per_gram_per_hour * np_saturation;

    // printf("An=%.2f  lim=%.2f  Nsat=%.2f  Psat=%.2f  adjPH=%.2f anConv=%f\n", An, limitation, input->nitrogen_saturation, input->phosphorus_saturation, adjusted_PH, input->convert_to_photosynthesis_per_gram_per_hour);
    // printf("adjusted: %f\n", adjusted_PH);
    return RMAX(REAL(0), RMIN(adjusted_PH, input->max_photosyntetic_rate));
    // return adjusted_PH > input->max_photosyntetic_rate ? input->max_photosyntetic_rate : adjusted_PH;
}

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
    real_t optimal_stechimetric_ratio,
    real_t phosphorus)
{
    real_t min_phosphorus_photosynthesis = min_nitrogen_photosynthesis / optimal_stechimetric_ratio;
    if (phosphorus < min_phosphorus_photosynthesis) {
        return 0;
    }
    return ((REAL(2.0) * phosphorus) / (phosphorus + min_phosphorus_photosynthesis)) - REAL(1.0);
}

real_t limitation_of_photosyntethic_rate(
    real_t starch,
    real_t feedback_on_photosynthesis,
    real_t max_starch)
{
    return feedback_on_photosynthesis + ((REAL(1.0) - feedback_on_photosynthesis) * ((max_starch - starch) / max_starch));
}

real_t max_starch(
    real_t max_starch_degradation_rate,
    real_t photoperiod)
{

    return max_starch_degradation_rate * (REAL(24.0) - photoperiod);
}

real_t starch_degradation(
    real_t max_starch_degradation_rate,
    real_t max_sucrose,
    real_t starch_night_start,
    real_t min_starch,
    real_t sucrose,
    real_t light,
    real_t time,
    real_t photoperiod)
{
    if (starch_night_start < min_starch) {
        return 0;
    }
    real_t degradation_rate = ((time / REAL(24.0)) + (REAL(1.0) - (time / REAL(24.0))) * (REAL(1.0) - (sucrose / (sucrose + max_sucrose)))) * ((starch_night_start - min_starch) / (REAL(24.0) - photoperiod));

    return (REAL(1.0) - light) * RMIN(max_starch_degradation_rate, degradation_rate);
}

real_t starch_sucrose_partition(
    real_t light,
    real_t starch_partition_coeff,
    real_t min_sucrose,
    real_t max_sucrose,
    real_t sucrose,
    real_t lambda_sdr,
    real_t lambda_sdi,
    real_t lambda_sni)
{
    // printf("\n light: %f, starch_part_coeff: %f, min_suc: %f, max_suc: %f, suc: %f, l_sdr: %f, l_sdi: %f, l_sni: %f\n", light,
    //  starch_partition_coeff,
    //  min_sucrose,
    //  max_sucrose,
    //  sucrose,
    //  lambda_sdr,
    //  lambda_sdi,
    //  lambda_sni);
    real_t s_min_fract = (min_sucrose / (min_sucrose + sucrose));
    real_t s_max_fract = (sucrose / (sucrose + max_sucrose));
    real_t day_part = light * ((-starch_partition_coeff * lambda_sdr * s_min_fract) + ((REAL(1.0) - starch_partition_coeff) * lambda_sdi * s_max_fract));

    real_t night_part = (REAL(1.0) - light) * (REAL(1.0) - starch_partition_coeff) * lambda_sni * s_min_fract;

    return day_part + night_part;
}
real_t starch_sucrose_partition_f(
    real_t t,
    real_t starch_partition_coeff,
    void* params)
{
    Input* input = (Input*)params;

    return starch_sucrose_partition(
        input->light,
        starch_partition_coeff,
        input->min_sucrose,
        input->max_sucrose,
        input->sucrose,
        input->lambda_sdr,
        input->lambda_sdi,
        input->lambda_sni);
}

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
    real_t min_nitrogen_photosynthesis)
{

    return nitrogen_uptake - ((respiration_frequency * (sucrose + starch) + (night_efficiency_starch * sucrose_loading_frequency)) * assimilaton_cost * (leaf_biomass / (total_biomass + epsylon))) - ((photosyntesis / max_photosyntetic_rate) * min_nitrogen_photosynthesis * nutrient_conversion_parameter);
}

real_t nitrogen_content_f(
    real_t t,
    real_t nitrogen,
    void* params)
{
    Input* input = (Input*)params;
    return nitrogen_content(
        input->nitrogen_uptake,
        input->respiration_frequency,
        input->sucrose,
        input->starch,
        input->sucrose_loading_frequency,
        input->night_efficiency_starch,
        input->assimilation_cost_nitrogen,
        input->leaf_biomass,
        input->total_biomass,
        input->photosynthesis,
        input->max_photosyntetic_rate,
        input->nutrient_conversion_parameter,
        input->min_nitrogen_photosynthesis);
}

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
    real_t min_phosphorus_photosynthesis)
{

    return phosphorus_uptake - ((respiration_frequency * (sucrose + starch) + (night_efficiency_starch * sucrose_loading_frequency)) * assimilaton_cost * (leaf_biomass / (total_biomass + epsylon))) - ((photosyntesis / max_photosyntetic_rate) * min_phosphorus_photosynthesis * nutrient_conversion_parameter);
}
real_t phosphorus_content_f(
    real_t t,
    real_t phosphorus,
    void* params)
{
    Input* input = (Input*)params;

    return phosphorus_content(
        input->phosphorus_uptake,
        input->respiration_frequency,
        input->sucrose,
        input->starch,
        input->sucrose_loading_frequency,
        input->night_efficiency_starch,
        input->assimilation_cost_phorphorus,
        input->leaf_biomass,
        input->total_biomass,
        input->photosynthesis,
        input->max_photosyntetic_rate,
        input->nutrient_conversion_parameter,
        input->min_phosphorus_photosynthesis);
}

real_t assimilation_cost_nitrogen(
    real_t lambda_csn)
{
    return lambda_csn * REAL(1.724);
}

real_t assimilation_cost_phosphorus(
    real_t assimilation_cost_nitrogen,
    real_t optimal_stechiometric_ratio)
{

    return assimilation_cost_nitrogen / optimal_stechiometric_ratio;
}

real_t pot_nitrogen_uptake(
    real_t nitrogen_nutrient_uptake,
    real_t nitrogen_affinity,
    real_t root_biomass,
    real_t total_biomass

)
{
    return nitrogen_nutrient_uptake * nitrogen_affinity * (root_biomass / (total_biomass + epsylon));
}

real_t nitrogen_uptake(
    real_t pot_nitrogen_uptake,
    real_t sucrose,
    real_t nitrogen_uptake_sucrose_consumption)
{

    return pot_nitrogen_uptake * (sucrose / (sucrose + (nitrogen_uptake_sucrose_consumption * pot_nitrogen_uptake * REAL(1.0)) + epsylon));
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
    real_t total_biomass

)
{
    return phosphorus_nutrient_uptake * phosphorus_affinity * (root_biomass / (total_biomass + epsylon));
}

real_t phosphorus_uptake(
    real_t pot_phosphorus_uptake,
    real_t sucrose,
    real_t phosphorus_uptake_sucrose_consumption)
{

    return pot_phosphorus_uptake * sucrose / (sucrose + (phosphorus_uptake_sucrose_consumption * pot_phosphorus_uptake * REAL(1.0)) + epsylon);
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
    real_t max_photosyntetic_rate,
    real_t lambda_k,
    real_t optimal_stochiometric_ratio,
    real_t phosphorus,
    real_t min_nitrogen,
    real_t nitrogen_uptake_sucrose_consumption,
    real_t starch_partition_coeff,
    real_t starch_degradation_rate)
{
    // printf(
    //  "\nnitrogen_affinity=%REAL(8.3)  nitrogen_uptake=%REAL(8.3)  nitrogen_cost=%REAL(8.3)  "
    //  "max_nitrogen=%REAL(8.3)  nitrogen=%REAL(8.3)  photosynthesis=%REAL(8.3)  "
    //  "max_photosyntetic_rate=%REAL(8.3)  lambda_k=%REAL(8.3)  optimal_stochiometric_ratio=%REAL(8.3)  "
    //  "phosphorus=%REAL(8.3)  min_nitrogen=%REAL(8.3)  nitrogen_uptake_sucrose_consumption=%REAL(8.3)  "
    //  "starch_partition_coeff=%REAL(8.3)  starch_degradation_rate=%REAL(8.3)\n",
    //  nitrogen_affinity,
    //  nitrogen_uptake,
    //  nitrogen_cost,
    //  max_nitrogen,
    //  nitrogen,
    //  photosynthesis,
    //  max_photosyntetic_rate,
    //  lambda_k,
    //  optimal_stochiometric_ratio,
    //  phosphorus,
    //  min_nitrogen,
    //  nitrogen_uptake_sucrose_consumption,
    //  starch_partition_coeff,
    //  starch_degradation_rate);

    real_t n_nmax = (max_nitrogen / (nitrogen + max_nitrogen));
    real_t un_real = (REAL(1.0) - (nitrogen_uptake / (nitrogen_uptake + nitrogen_cost + epsylon)));
    real_t n_op = nitrogen_affinity * lambda_k * ((nitrogen - (optimal_stochiometric_ratio * phosphorus)) / (nitrogen + (optimal_stochiometric_ratio * phosphorus) + epsylon));

    real_t term1 = (REAL(1.0) - nitrogen_affinity) * ((un_real * n_nmax) + (photosynthesis / max_photosyntetic_rate) - n_op);

    real_t n_nmin = (nitrogen / (nitrogen + min_nitrogen));
    real_t term2 = nitrogen_affinity * (n_nmin + ((nitrogen_uptake_sucrose_consumption * nitrogen_uptake) / ((nitrogen_uptake_sucrose_consumption * nitrogen_uptake) + (photosynthesis * (REAL(1.0) - starch_partition_coeff)) + starch_degradation_rate + epsylon)));
    // //printf("\nn_nmax: %f, un_real %f, n_op %f, term1: %f, n_nmin: %f, term2: %f, term: %f\n", n_nmax, un_real, n_op, term1, n_nmin, term2, term1 - term2);
    return term1 - term2;
}

real_t nitrogen_affinity_f(
    real_t t,
    real_t nitrogen_affinity_p,
    void* params)
{
    Input* in = (Input*)params;
    // real_t loc_nitrogen_nutrient_uptake = nitrogen_nutrient_uptake(input->max_nitrogen_uptake, input->nitrogen_soil_content, input->Michaelis_Menten_constant_nitrogen);
    // real_t loc_pot_nitrogen_uptake = pot_nitrogen_uptake(loc_nitrogen_nutrient_uptake, nitrogen_affinity_p, input->root_biomass, input->total_biomass);
    // real_t loc_nitrogen_uptake = nitrogen_uptake(loc_pot_nitrogen_uptake, input->sucrose, input->nitrogen_uptake_sucrose_consumption);
    // real_t loc_min_nitrogen = min_nitrogen(input->respiration_frequency, input->max_sucrose, input->max_starch, input->sucrose_loading_frequency, input->assimilation_cost_nitrogen, input->min_nitrogen_photosynthesis, input->nutrient_conversion_parameter, input->photoperiod);
    // real_t loc_max_nitrogen = max_nitrogen(loc_min_nitrogen);
    // real_t loc_nitrogen_cost = nitrogen_cost(input->respiration_frequency, input->sucrose, input->starch, input->sucrose_loading_frequency, input->night_efficiency_starch, input->assimilation_cost_nitrogen, input->leaf_biomass, input->total_biomass, input->photosynthesis, input->max_photosyntetic_rate, input->min_nitrogen_photosynthesis);

    return nitrogen_affinity(nitrogen_affinity_p,
        in->nitrogen_uptake,
        in->nitrogen_cost,
        in->max_nitrogen,
        in->nitrogen,
        in->photosynthesis,
        in->max_photosyntetic_rate,
        in->lambda_k,
        in->optimal_stechiometric_ratio,
        in->phosphorus,
        in->min_nitrogen,
        in->nitrogen_uptake_sucrose_consumption,
        in->starch_partition_coeff,
        in->starch_degradation_rate);
}

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
    real_t starch_degradation_rate)
{
    // //printf("\np_aff: %f ", phosphorus_affinity);
    // printf(
    // "\nphosphorus_affinity=%f "
    // "phosphorus_uptake=%f "
    // "phosphorus_cost=%f "
    // "max_phosphorus=%f "
    // "phosphorus=%f "
    // "photosynthesis=%f "
    // "max_photosyntetic_rate=%f "
    // "lambda_k=%f "
    // "optimal_stochiometric_ratio=%f "
    // "nitrogen=%f "
    // "min_phosphorus=%f "
    // "phosphorus_uptake_sucrose_consumption=%f "
    // "starch_partition_coeff=%f "
    // "starch_degradation_rate=%f\n",
    // phosphorus_affinity,
    // phosphorus_uptake,
    // phosphorus_cost,
    // max_phosphorus,
    // phosphorus,
    // photosynthesis,
    // max_photosyntetic_rate,
    // lambda_k,
    // optimal_stochiometric_ratio,
    // nitrogen,
    // min_phosphorus,
    // phosphorus_uptake_sucrose_consumption,
    // starch_partition_coeff,
    // starch_degradation_rate);

    real_t n_nmax = (max_phosphorus / (phosphorus + max_phosphorus));
    real_t un_real = (REAL(1.0) - (phosphorus_uptake / (phosphorus_uptake + phosphorus_cost + epsylon)));
    real_t n_op = phosphorus_affinity * lambda_k * ((nitrogen - (optimal_stochiometric_ratio * phosphorus)) / (nitrogen + (optimal_stochiometric_ratio * phosphorus) + epsylon));

    real_t term1 = (REAL(1.0) - phosphorus_affinity) * ((un_real * n_nmax) + (photosynthesis / max_photosyntetic_rate) + n_op);

    real_t n_nmin = (phosphorus / (phosphorus + min_phosphorus));
    real_t term2 = phosphorus_affinity * (n_nmin + ((phosphorus_uptake_sucrose_consumption * phosphorus_uptake) / ((phosphorus_uptake_sucrose_consumption * phosphorus_uptake) + (photosynthesis * (REAL(1.0) - starch_partition_coeff)) + starch_degradation_rate + epsylon)));
    // //printf("\np_nmax: %f, up_real %f, p_op %f, term1: %f, p_nmin: %f, term2: %f, term: %f\n", n_nmax, un_real, n_op, term1, n_nmin, term2, term1 - term2);
    return term1 - term2;
}
real_t phosphorus_affinity_f(
    real_t t,
    real_t phosphorus_affinity_p,
    void* params)
{
    Input* input = (Input*)params;
    // real_t loc_phosphorus_nutrient_uptake = phosphorus_nutrient_uptake(input->max_phosphorus_uptake, input->phosphorus_soil_content, input->Michaelis_Menten_constant_phosphorus);
    // real_t loc_pot_phosphorus_uptake = pot_phosphorus_uptake(loc_phosphorus_nutrient_uptake, phosphorus_affinity_p, input->root_biomass, input->total_biomass);
    // real_t loc_phosphorus_uptake = phosphorus_uptake(loc_pot_phosphorus_uptake, input->sucrose, input->phosphorus_uptake_sucrose_consumption);
    //
    // real_t loc_min_phosphorus = min_phosphorus(
    //     input->respiration_frequency, input->max_sucrose, input->max_starch, input->sucrose_loading_frequency, input->assimilation_cost_phorphorus, input->min_phosphorus_photosynthesis, input->nutrient_conversion_parameter, input->photoperiod);
    //
    // real_t loc_max_phosphorus = max_phosphorus(input->max_nitrogen, input->optimal_stechiometric_ratio);

    return phosphorus_affinity(phosphorus_affinity_p, input->phosphorus_uptake, input->phosphorus_cost, input->max_phosphorus, input->phosphorus, input->photosynthesis, input->max_photosyntetic_rate, input->lambda_k, input->optimal_stechiometric_ratio, input->nitrogen, input->min_phosphorus, input->phosphorus_uptake_sucrose_consumption, input->starch_partition_coeff, input->starch_degradation_rate);
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
    real_t max_photosyntetic_rate,
    real_t min_nitrogen_photosyntesis)
{
    real_t first_term = (respiration_frequency * (sucrose + starch) + (night_efficiency_starch * sucrose_loading_frequency))
        * assimilation_cost_nitrogen * (leaf_biomass / (total_biomass + epsylon));

    real_t second_term = (photosynthesis / max_photosyntetic_rate) * min_nitrogen_photosyntesis;

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
    real_t max_photosyntetic_rate,
    real_t min_phosphorus_photosyntesis)
{
    real_t first_term = (respiration_frequency * (sucrose + starch) + (night_efficiency_starch * sucrose_loading_frequency))
        * assimilation_cost_phosphorus * (leaf_biomass / (total_biomass + epsylon));

    real_t second_term = (photosynthesis / max_photosyntetic_rate) * min_phosphorus_photosyntesis;

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

real_t max_phosphorus(real_t max_nitrogen, real_t optimal_stechiometric_ratio)
{
    return max_nitrogen / optimal_stechiometric_ratio;
}

real_t starch_production(
    real_t photosynthesis,
    real_t starch_partition_coeff,
    real_t starch_degradation_rate)
{

    return (starch_partition_coeff * photosynthesis) - starch_degradation_rate;
}
real_t starch_production_f(
    real_t t,
    real_t starch,
    void* params)
{
    Input* input = (Input*)params;
    return starch_production(input->photosynthesis, input->starch_partition_coeff, input->starch_degradation_rate);
}

real_t sucrose_production(
    real_t sucrose_partition_coeff,
    real_t photosynthesis,
    real_t starch_degradation_rate,
    real_t uptake_cost,
    real_t transport_cost,
    real_t respiration_frequency,
    real_t sucrose,
    real_t sucrose_loading_frequency,
    real_t night_efficiency_starch)
{
    // printf(
    //  "sucrose_partition_coeff=%REAL(8.3)  photosynthesis=%REAL(8.3)  starch_degradation_rate=%REAL(8.3)  "
    //  "uptake_cost=%REAL(8.3)  transport_cost=%REAL(8.3)  respiration_frequency=%REAL(8.3)  "
    //  "sucrose=%REAL(8.3)  sucrose_loading_frequency=%REAL(8.3)  night_efficiency_starch=%REAL(8.3)\n",
    //  sucrose_partition_coeff,
    //  photosynthesis,
    //  starch_degradation_rate,
    //  uptake_cost,
    //  transport_cost,
    //  respiration_frequency,
    //  sucrose,
    //  sucrose_loading_frequency,
    //  night_efficiency_starch);

    return ((REAL(1.0) - sucrose_partition_coeff) * photosynthesis) + starch_degradation_rate - uptake_cost - transport_cost - (respiration_frequency * sucrose) - (sucrose_loading_frequency * night_efficiency_starch);
}
real_t sucrose_production_f(
    real_t t,
    real_t sucrose,
    void* params)
{
    Input* in = (Input*)params;

    // real_t loc_photosynthesis = photosynthesis(
    //     in->light,
    //     in->limitation_of_photosyntetic_rate,
    //     in->max_photosyntetic_rate,
    //     in->nitrogen_saturation,
    //     in->phosphorus_saturation,
    //     in->leaf_biomass,
    //     in->min_leaf_biomass);
    // real_t loc_night_efficiency_starch = night_efficieny_starch(sucrose, in->max_sucrose, in->lambda_g, in->light, in->starch_partition_coeff);
    // real_t loc_uptake_cost = uptake_cost(in->nitrogen_uptake_sucrose_consumption,
    //     in->nitrogen_uptake,
    //     in->phosphorus_uptake_sucrose_consumption,
    //     in->phosphorus_uptake);
    //
    // real_t loc_transport_cost = transport_cost(in->sucrose_consumption_transport,
    //     in->respiration_frequency,
    //     sucrose,
    //     in->sucrose_loading_frequency,
    //     loc_night_efficiency_starch);
    //
    // return sucrose_production(in->starch_partition_coeff,
    //     loc_photosynthesis,
    //     in->starch_degradation_rate,
    //     loc_uptake_cost,
    //     loc_transport_cost,
    //     in->respiration_frequency,
    //     sucrose,
    //     in->sucrose_loading_frequency,
    //     loc_night_efficiency_starch);
    return sucrose_production(in->starch_partition_coeff,
        in->photosynthesis,
        in->starch_degradation_rate,
        in->uptake_cost,
        in->transport_cost,
        in->respiration_frequency,
        sucrose,
        in->sucrose_loading_frequency,
        in->night_efficiency_starch);
}

real_t night_efficieny_starch(
    real_t sucrose,
    real_t max_sucrose,
    real_t lambda_g,
    real_t light,
    real_t starch_partition_coeff)
{

    return RMIN(sucrose, max_sucrose) * (lambda_g + ((REAL(1.0) - lambda_g) * (REAL(1.0) - light + (light * (REAL(1.0) - starch_partition_coeff)))));
}

real_t uptake_cost(
    real_t nitrogen_uptake_sucrose_consumption,
    real_t nitrogen_uptake,
    real_t phosphorus_uptake_sucrose_consumption,
    real_t phosphorus_uptake)
{
    return (nitrogen_uptake_sucrose_consumption * nitrogen_uptake) + (phosphorus_uptake_sucrose_consumption * phosphorus_uptake);
}

real_t transport_cost(
    real_t sucrose_consumption_transport,
    real_t respiration_frequency,
    real_t sucrose,
    real_t sucrose_loading_frequency,
    real_t night_efficiency_starch)
{
    return sucrose_consumption_transport * ((respiration_frequency * sucrose) + (sucrose_loading_frequency * night_efficiency_starch));
}

real_t leaf_growth(
    real_t labmda_sb,
    real_t sucrose_root_allocation,
    real_t sucrose_loading_frequency,
    real_t night_efficiency_starch,
    real_t leaf_biomass,
    real_t leaf_death_rate,
    real_t leaf_competative_rate)
{
    // printf("\n[LEAF GROWTH DEBUG] lambda_sb: %.6f", labmda_sb);
    // printf(" sucrose_root_allocation: %.6f", sucrose_root_allocation);
    // printf(" sucrose_loading_frequency: %.6f", sucrose_loading_frequency);
    // printf(" night_efficiency_starch: %.6f", night_efficiency_starch);
    // printf(" leaf_biomass: %.6f", leaf_biomass);
    // printf(" leaf_death_rate: %.6f", leaf_death_rate);
    // printf(" leaf_competative_rate: %.6f", leaf_competative_rate);

    return (labmda_sb * (REAL(1.0) - sucrose_root_allocation) * sucrose_loading_frequency * night_efficiency_starch * leaf_biomass) - (leaf_death_rate * leaf_biomass) - (leaf_competative_rate * leaf_biomass * leaf_biomass);
}
real_t leaf_growth_f(
    real_t t,
    real_t leaf_biomass,
    void* params)
{
    Input* input = (Input*)params;

    return leaf_growth(
        input->lambda_sb,
        input->sucrose_root_allocation,
        input->sucrose_loading_frequency,
        input->night_efficiency_starch,
        leaf_biomass,
        input->leaf_deathrate,
        input->leaf_competitive_rate);
}

real_t root_growth(
    real_t labmda_sb,
    real_t sucrose_root_allocation,
    real_t sucrose_loading_frequency,
    real_t night_efficiency_starch,
    real_t leaf_biomass,
    real_t root_biomass,
    real_t root_death_rate,
    real_t root_competative_rate)
{

    // printf(
    //     "\n[ROOT GROWTH DEBUG]"
    //     " lambda_sb= %.6f"
    //     " sucrose_root_allocation= %.6f"
    //     " sucrose_loading_frequency= %.6f"
    //     " night_efficiency_starch= %.6f"
    //     " leaf_biomass= %.6f"
    //     " root_biomass= %.6f"
    //     " root_death_rate= %.6f"
    //     " root_competative_rate= %.6f",
    //     labmda_sb,
    //     sucrose_root_allocation,
    //     sucrose_loading_frequency,
    //     night_efficiency_starch,
    //     leaf_biomass,
    //     root_biomass,
    //     root_death_rate,
    //     root_competative_rate);
    return (labmda_sb * sucrose_root_allocation * sucrose_loading_frequency * night_efficiency_starch * leaf_biomass) - (root_death_rate * root_biomass) - (root_competative_rate * root_biomass * root_biomass);
}
real_t root_growth_f(
    real_t t,
    real_t root_biomass,
    void* params)
{

    Input* input = (Input*)params;
    return root_growth(
        input->lambda_sb,
        input->sucrose_root_allocation,
        input->sucrose_loading_frequency,
        input->night_efficiency_starch,
        input->leaf_biomass,
        root_biomass,
        input->root_deathrate,
        input->root_competitive_rate);
}

real_t stochiometric_signal(
    real_t optimal_stechiometric_ratio,
    real_t nitrogen,
    real_t phosphorus)
{
    return optimal_stechiometric_ratio / (optimal_stechiometric_ratio + (nitrogen / (phosphorus + epsylon)));
}

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
    real_t min_phosphorus)
{
    // printf(
    // "\n[SUCROSE ROOT ALLOCATION] sucrose_root_allocation=%REAL(8.3)  nitrogen_affinity=%REAL(8.3)  phospsorus_affinity=%REAL(8.3)  "
    // "stochiometric_signal=%REAL(8.3)  sucrose=%REAL(8.3)  min_sucrose=%REAL(8.3)  "
    // "nitrogen=%REAL(8.3)  min_nitrogen=%REAL(8.3)  phosphorus=%REAL(8.3)  min_phosphorus=%REAL(8.3)\n",
    // sucrose_root_allocation,
    // nitrogen_affinity,
    // phospsorus_affinity,
    // stochiometric_signal,
    // sucrose,
    // min_sucrose,
    // nitrogen,
    // min_nitrogen,
    // phosphorus,
    // min_phosphorus);

    real_t first = (REAL(1.0) - sucrose_root_allocation) * ((nitrogen_affinity * stochiometric_signal) + ((REAL(1.0) - stochiometric_signal) * phospsorus_affinity));
    real_t secound = sucrose_root_allocation * (((nitrogen * stochiometric_signal) / (nitrogen + min_nitrogen)) + ((phosphorus * (REAL(1.0) - stochiometric_signal)) / (phosphorus + min_phosphorus)) + (min_sucrose / (min_sucrose + sucrose)));
    return first - secound;
}
real_t sucrose_root_allocation_f(
    real_t t,
    real_t sucrose_root_allocation_p,
    void* params)
{
    Input* input = (Input*)params;
    return sucrose_root_allocation(
        sucrose_root_allocation_p,
        input->nitrogen_affinity,
        input->phosphorus_affinity,
        input->stochiometric_signal,
        input->sucrose,
        input->min_sucrose,
        input->nitrogen,
        input->min_nitrogen,
        input->phosphorus,
        input->min_phosphorus);
}

// enviromental contrubution tedone

real_t lambda_sb_f(real_t lambda_sb, real_t nitrogen_soil, real_t phosphorus_soil)
{
    real_t n = nitrogen_soil / REAL(12.0);
    real_t p = phosphorus_soil / REAL(0.15);
    real_t g = -(REAL(0.5) * (RPOW(n, 2) + RPOW(p, 2))) + (n + p);
    real_t R_minus = RMAX(g, 0);
    if (RMIN(n, p) >= 1) {
        return lambda_sb;
    } else {
        return lambda_sb * R_minus;
    }
}

real_t R_plus(real_t nitrogen_soil, real_t phosphorus_soil)
{
    real_t n = nitrogen_soil / REAL(12.0);
    real_t p = phosphorus_soil / REAL(0.15);

    real_t R_plus = (REAL(69.0) * (RPOW(n, 2) + RPOW(p, 2))) - REAL(138.0) * (n + p) + REAL(139);
    if (RMIN(n, p) >= 1) {
        return 1;
    }
    return R_plus;
}

// Farquhar C3 Photosynthesis Model

real_t rubisco_limited_photosynthesis(
    real_t intercellular_CO2, // Ci
    real_t CO2_compensation_point, // Gamma_star
    real_t max_rubisco_carboxylation_rate, // Vcmax
    real_t michaelis_constant_CO2, // Kc
    real_t michaelis_constant_O2, // Ko
    real_t oxygen_concentration // O
)
{
    real_t denominator = intercellular_CO2 + (michaelis_constant_CO2 * (REAL(1.0) + (oxygen_concentration / michaelis_constant_O2)));
    return max_rubisco_carboxylation_rate * (intercellular_CO2 - CO2_compensation_point) / denominator;
}

// Light-limited (RuBP regeneration) assimilation rate
real_t light_limited_photosynthesis(
    real_t intercellular_CO2, // Ci
    real_t CO2_compensation_point, // Gamma_star
    real_t electron_transport_rate // J
)
{
    real_t denominator = (REAL(4.5) * intercellular_CO2) + (REAL(10.5) * CO2_compensation_point);
    return electron_transport_rate * (intercellular_CO2 - CO2_compensation_point) / denominator;
}

// Net photosynthesis rate
real_t net_photosynthesis(
    real_t rubisco_limited_rate, // Ac
    real_t light_limited_rate, // Aj
    real_t respiration_rate // Rd
)
{
    real_t limiting_rate = RMIN(rubisco_limited_rate, light_limited_rate);
    return limiting_rate - respiration_rate;
}

real_t convert_to_photosynthesis_per_gram_per_hour(
    real_t net_photosynthesis_rate, // An
    real_t leaf_biomass,
    real_t specific_leaf_area)
{

    return net_photosynthesis_rate * REAL(0.0025) * 600;
}

// Intercellular CO2 based on diffusion from ambient and net photosynthesis
real_t intercellular_CO2(
    real_t ambient_CO2_concentration, // Ca
    real_t net_photosynthesis_rate, // An
    real_t stomatal_conductance // gs
)
{
    return ambient_CO2_concentration - (net_photosynthesis_rate / stomatal_conductance);
}

real_t iterate_ci(Input* input, int max_iter, real_t epsilon)
{
    real_t Ci = input->ambient_CO2_concentration * REAL(0.7); // inicijalna pretpostavka
    real_t prev_An = 0.0;
    real_t An = 0.0;

    for (int i = 0; i < max_iter; i++) {
        real_t Ac = rubisco_limited_photosynthesis(
            Ci,
            input->CO2_compensation_point,
            input->max_rubisco_carboxylation_rate,
            input->michaelis_constant_CO2,
            input->michaelis_constant_O2,
            input->oxygen_concentration);

        real_t Aj = light_limited_photosynthesis(
            Ci,
            input->CO2_compensation_point,
            input->electron_transport_rate);

        An = net_photosynthesis(Ac, Aj, input->respiration_rate);

        real_t new_Ci = input->ambient_CO2_concentration - An / input->stomatal_conductance;

        if (RABS(new_Ci - Ci) < epsilon && RABS(An - prev_An) < epsilon) {
            break;
        }

        Ci = new_Ci;
        prev_An = An;
    }

    input->intercellular_CO2 = Ci;
    input->net_photosynthesis = An;
    input->net_photosynthesis_rate = An;
    return An;
}

// Jarvis stomal

real_t calculate_stomatal_conductance_jarvis(
    real_t incoming_PAR,
    real_t vapor_pressure_deficit,
    real_t leaf_temperature_celsius,
    real_t ambient_CO2_concentration,
    real_t leaf_water_potential)
{
    const real_t b1 = REAL(0.6);
    const real_t b2 = REAL(0.002);
    const real_t b10 = REAL(0.05);
    real_t q = b10 / b1;
    real_t f_PAR = (incoming_PAR <= q) ? b10 : (b1 * b2 * (incoming_PAR - q)) / (b1 + b2 * (incoming_PAR - q));

    const real_t b5 = REAL(0.2);
    real_t f_VPD = RMAX(REAL(0.0), REAL(1.0) - b5 * vapor_pressure_deficit);

    const real_t T0 = REAL(25.0);
    const real_t T1 = REAL(5.0);
    const real_t Th = REAL(45.0);
    const real_t b4 = (Th - T0) / (Th - T1);
    const real_t b3 = REAL(1.0) / ((T0 - T1) * RPOW((Th - T1), b4));
    real_t f_T = b3 * (leaf_temperature_celsius - T1) * RPOW((Th - leaf_temperature_celsius), b4);
    f_T = RMAX(REAL(0.0), RMIN(f_T, REAL(1.0)));

    const real_t b8 = REAL(0.2);
    real_t f_CO2 = REAL(1.0);
    if (ambient_CO2_concentration < 100.0) {
        f_CO2 = 1.0;
    } else if (ambient_CO2_concentration < 1000.0) {
        const real_t b7 = (1.0 - b8) / 900.0;
        f_CO2 = 1.0 - b7 * (ambient_CO2_concentration - 100.0);
    } else {
        f_CO2 = b8;
    }

    const real_t b6 = 3.0;
    const real_t psi_min = -1.5;
    real_t sigma_psi = leaf_water_potential - psi_min;
    real_t f_water = 1.0 - REXP(-b6 * sigma_psi);
    f_water = RMAX(0.0, RMIN(f_water, 1.0));

    return f_PAR * f_VPD * f_T * f_CO2 * f_water;
}

real_t calculate_vapor_pressure_deficit(real_t temperature_celsius, real_t relative_humidity_percent)
{
    real_t saturation_vapor_pressure = REAL(0.61078) * REXP((REAL(17.27) * temperature_celsius) / (temperature_celsius + REAL(237.3)));
    real_t actual_vapor_pressure = saturation_vapor_pressure * (relative_humidity_percent / REAL(100.0));
    return saturation_vapor_pressure - actual_vapor_pressure;
}

void update_light_conditions(Input* input, real_t current_hour)
{
    const real_t max_PAR = 1200.0;

    real_t sunrise = 0.0;
    real_t sunset = sunrise + input->photoperiod;

    if (current_hour >= sunrise && current_hour < sunset) {
        input->light = 1.0;
        real_t day_fraction = (current_hour - sunrise) / (sunset - sunrise);
        input->light_PAR = max_PAR * RSIN(day_fraction * M_PI);
    } else {
        input->light = 0.0;
        input->light_PAR = 0.0;
    }
}

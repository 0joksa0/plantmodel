#include "model/model.h"
#include "simulation/solver/solver.h"
#include <math.h>
#include <stdio.h>

float epsylon = 0.000001f;
float photosynthesis(
    float light,
    float limitation_of_photosynthetic_rate,
    float max_photosyntetic_rate,
    float nitrogen_saturation,
    float phosphorus_saturation,
    float leaf_biomass,
    float min_leaf_biomass)
{
    if (leaf_biomass < min_leaf_biomass) {
        return 0;
    }

    //printf("lim_p_rate: %f, nit_sat: %f, phos_sat: %f\n", limitation_of_photosynthetic_rate, nitrogen_saturation, phosphorus_saturation);
    float photo = light * limitation_of_photosynthetic_rate * max_photosyntetic_rate * fminf(nitrogen_saturation, phosphorus_saturation);
    return photo;
}

float nitrogen_saturation(
    float nitrogen,
    float min_nitrogen_photosynthesis)
{
    if (nitrogen < min_nitrogen_photosynthesis) {
        return 0;
    }
    return ((2.0f * nitrogen) / (nitrogen + min_nitrogen_photosynthesis)) - 1.0f;
}

float phosphorus_saturation(
    float min_nitrogen_photosynthesis,
    float optimal_stechimetric_ratio,
    float phosphorus)
{
    float min_phosphorus_photosynthesis = min_nitrogen_photosynthesis / optimal_stechimetric_ratio;
    if (phosphorus < min_phosphorus_photosynthesis) {
        return 0;
    }
    return ((2.0f * phosphorus) / (phosphorus + min_phosphorus_photosynthesis)) - 1.0f;
}

float limitation_of_photosyntethic_rate(
    float starch,
    float feedback_on_photosynthesis,
    float max_starch)
{
    return feedback_on_photosynthesis + ((1.0f - feedback_on_photosynthesis) * ((max_starch - starch) / max_starch));
}

float max_starch(
    float max_starch_degradation_rate,
    float photoperiod)
{

    return max_starch_degradation_rate * (24.0f - photoperiod);
}

float starch_degradation(
    float max_starch_degradation_rate,
    float max_sucrose,
    float starch_night_start,
    float min_starch,
    float sucrose,
    float light,
    float time,
    float photoperiod)
{
    if (starch_night_start < min_starch) {
        return 0;
    }
    float degradation_rate = ((time / 24.0f) + (1.0f - (time / 24.0f)) * (1.0f - (sucrose / (sucrose + max_sucrose)))) * ((starch_night_start - min_starch) / (24.0f - photoperiod));

    return (1.0f - light) * fminf(max_starch_degradation_rate, degradation_rate);
}

float starch_sucrose_partition(
    float light,
    float starch_partition_coeff,
    float min_sucrose,
    float max_sucrose,
    float sucrose,
    float lambda_sdr,
    float lambda_sdi,
    float lambda_sni)
{
    //printf("\n light: %f, starch_part_coeff: %f, min_suc: %f, max_suc: %f, suc: %f, l_sdr: %f, l_sdi: %f, l_sni: %f\n", light,
        // starch_partition_coeff,
        // min_sucrose,
        // max_sucrose,
        // sucrose,
        // lambda_sdr,
        // lambda_sdi,
        // lambda_sni);
    float s_min_fract = (min_sucrose / (min_sucrose + sucrose));
    float s_max_fract = (sucrose / (sucrose + max_sucrose));
    float day_part = light * ((-starch_partition_coeff * lambda_sdr * s_min_fract) + ((1.0f - starch_partition_coeff) * lambda_sdi * s_max_fract));

    float night_part = (1.0f - light) * (1.0f - starch_partition_coeff) * lambda_sni * s_min_fract;

    return day_part + night_part;
}
float starch_sucrose_partition_f(
    float t,
    float starch_partition_coeff,
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
    float min_nitrogen_photosynthesis)
{

    return nitrogen_uptake - ((respiration_frequency * (sucrose + starch) + (night_efficiency_starch * sucrose_loading_frequency)) * assimilaton_cost * (leaf_biomass / (total_biomass + epsylon))) - ((photosyntesis / max_photosyntetic_rate) * min_nitrogen_photosynthesis * nutrient_conversion_parameter);
}

float nitrogen_content_f(
    float t,
    float nitrogen,
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
    float min_phosphorus_photosynthesis)
{

    return phosphorus_uptake - ((respiration_frequency * (sucrose + starch) + (night_efficiency_starch * sucrose_loading_frequency)) * assimilaton_cost * (leaf_biomass / (total_biomass + epsylon))) - ((photosyntesis / max_photosyntetic_rate) * min_phosphorus_photosynthesis * nutrient_conversion_parameter);
}
float phosphorus_content_f(
    float t,
    float phosphorus,
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

float assimilation_cost_nitrogen(
    float lambda_csn)
{
    return lambda_csn * 1.724f;
}

float assimilation_cost_phosphorus(
    float assimilation_cost_nitrogen,
    float optimal_stechiometric_ratio)
{

    return assimilation_cost_nitrogen / optimal_stechiometric_ratio;
}

float pot_nitrogen_uptake(
    float nitrogen_nutrient_uptake,
    float nitrogen_affinity,
    float root_biomass,
    float total_biomass

)
{
    return nitrogen_nutrient_uptake * nitrogen_affinity * (root_biomass / (total_biomass + epsylon));
}

float nitrogen_uptake(
    float pot_nitrogen_uptake,
    float sucrose,
    float nitrogen_uptake_sucrose_consumption)
{

    return pot_nitrogen_uptake * (sucrose / (sucrose + (nitrogen_uptake_sucrose_consumption * pot_nitrogen_uptake * 1.0f) + epsylon));
}

float nitrogen_nutrient_uptake(
    float max_nitrogen_uptake,
    float nitrogen_soil_content,
    float Michaelis_Menten_const_nitrogen)
{
    return max_nitrogen_uptake * (nitrogen_soil_content / (nitrogen_soil_content + Michaelis_Menten_const_nitrogen));
}

float pot_phosphorus_uptake(
    float phosphorus_nutrient_uptake,
    float phosphorus_affinity,
    float root_biomass,
    float total_biomass

)
{
    return phosphorus_nutrient_uptake * phosphorus_affinity * (root_biomass / (total_biomass + epsylon));
}

float phosphorus_uptake(
    float pot_phosphorus_uptake,
    float sucrose,
    float phosphorus_uptake_sucrose_consumption)
{

    return pot_phosphorus_uptake * sucrose / (sucrose + (phosphorus_uptake_sucrose_consumption * pot_phosphorus_uptake * 1.0f) + epsylon);
}

float phosphorus_nutrient_uptake(
    float max_phosphorus_uptake,
    float phosphorus_soil_content,
    float Michaelis_Menten_const_phosphorus)
{
    return max_phosphorus_uptake * (phosphorus_soil_content / (phosphorus_soil_content + Michaelis_Menten_const_phosphorus));
}
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
    float starch_degradation_rate)
{
    //printf(
        // "\nnitrogen_affinity=%8.3f  nitrogen_uptake=%8.3f  nitrogen_cost=%8.3f  "
        // "max_nitrogen=%8.3f  nitrogen=%8.3f  photosynthesis=%8.3f  "
        // "max_photosyntetic_rate=%8.3f  lambda_k=%8.3f  optimal_stochiometric_ratio=%8.3f  "
        // "phosphorus=%8.3f  min_nitrogen=%8.3f  nitrogen_uptake_sucrose_consumption=%8.3f  "
        // "starch_partition_coeff=%8.3f  starch_degradation_rate=%8.3f\n",
        // nitrogen_affinity,
        // nitrogen_uptake,
        // nitrogen_cost,
        // max_nitrogen,
        // nitrogen,
        // photosynthesis,
        // max_photosyntetic_rate,
        // lambda_k,
        // optimal_stochiometric_ratio,
        // phosphorus,
        // min_nitrogen,
        // nitrogen_uptake_sucrose_consumption,
        // starch_partition_coeff,
        // starch_degradation_rate);

    float n_nmax = (max_nitrogen / (nitrogen + max_nitrogen));
    float un_real = (1.0f - (nitrogen_uptake / (nitrogen_uptake + nitrogen_cost + epsylon)));
    float n_op = nitrogen_affinity * lambda_k * ((nitrogen - (optimal_stochiometric_ratio * phosphorus)) / (nitrogen + (optimal_stochiometric_ratio * phosphorus) + epsylon));

    float term1 = (1.0f - nitrogen_affinity) * ((un_real * n_nmax) + (photosynthesis / max_photosyntetic_rate) - n_op);

    float n_nmin = (nitrogen / (nitrogen + min_nitrogen));
    float term2 = nitrogen_affinity * (n_nmin + ((nitrogen_uptake_sucrose_consumption * nitrogen_uptake) / ((nitrogen_uptake_sucrose_consumption * nitrogen_uptake) + (photosynthesis * (1.0f - starch_partition_coeff)) + starch_degradation_rate + epsylon)));
    // //printf("\nn_nmax: %f, un_real %f, n_op %f, term1: %f, n_nmin: %f, term2: %f, term: %f\n", n_nmax, un_real, n_op, term1, n_nmin, term2, term1 - term2);
    return term1 - term2;
}

float nitrogen_affinity_f(
    float t,
    float nitrogen_affinity_p,
    void* params)
{
    Input* input = (Input*)params;
    float loc_nitrogen_nutrient_uptake = nitrogen_nutrient_uptake(input->max_nitrogen_uptake, input->nitrogen_soil_content, input->Michaelis_Menten_constant_nitrogen);
    float loc_pot_nitrogen_uptake = pot_nitrogen_uptake(loc_nitrogen_nutrient_uptake, nitrogen_affinity_p, input->root_biomass, input->total_biomass);
    float loc_nitrogen_uptake = nitrogen_uptake(loc_pot_nitrogen_uptake, input->sucrose, input->nitrogen_uptake_sucrose_consumption);
    float loc_min_nitrogen = min_nitrogen(input->respiration_frequency, input->max_sucrose, input->max_starch, input->sucrose_loading_frequency, input->assimilation_cost_nitrogen, input->min_nitrogen_photosynthesis, input->nutrient_conversion_parameter, input->photoperiod);
    float loc_max_nitrogen = max_nitrogen(loc_min_nitrogen);
    float loc_nitrogen_cost = nitrogen_cost(input->respiration_frequency, input->sucrose, input->starch, input->sucrose_loading_frequency, input->night_efficiency_starch, input->assimilation_cost_nitrogen, input->leaf_biomass, input->total_biomass, input->photosynthesis, input->max_photosyntetic_rate, input->min_nitrogen_photosynthesis);

    return nitrogen_affinity(nitrogen_affinity_p,
        loc_nitrogen_uptake,
        loc_nitrogen_cost,
        loc_max_nitrogen,
        input->nitrogen,
        input->photosynthesis,
        input->max_photosyntetic_rate,
        input->lambda_k,
        input->optimal_stechiometric_ratio,
        input->phosphorus,
        loc_min_nitrogen,
        input->nitrogen_uptake_sucrose_consumption,
        input->starch_partition_coeff,
        input->starch_degradation_rate);
}

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
    float starch_degradation_rate)
{
    // //printf("\np_aff: %f ", phosphorus_affinity);
    //printf(
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

    float n_nmax = (max_phosphorus / (phosphorus + max_phosphorus));
    float un_real = (1.0f - (phosphorus_uptake / (phosphorus_uptake + phosphorus_cost + epsylon)));
    float n_op = phosphorus_affinity * lambda_k * ((nitrogen - (optimal_stochiometric_ratio * phosphorus)) / (nitrogen + (optimal_stochiometric_ratio * phosphorus) + epsylon));

    float term1 = (1.0f - phosphorus_affinity) * ((un_real * n_nmax) + (photosynthesis / max_photosyntetic_rate) + n_op);

    float n_nmin = (phosphorus / (phosphorus + min_phosphorus));
    float term2 = phosphorus_affinity * (n_nmin + ((phosphorus_uptake_sucrose_consumption * phosphorus_uptake) / ((phosphorus_uptake_sucrose_consumption * phosphorus_uptake) + (photosynthesis * (1.0f - starch_partition_coeff)) + starch_degradation_rate + epsylon)));
    // //printf("\np_nmax: %f, up_real %f, p_op %f, term1: %f, p_nmin: %f, term2: %f, term: %f\n", n_nmax, un_real, n_op, term1, n_nmin, term2, term1 - term2);
    return term1 - term2;
}
float phosphorus_affinity_f(
    float t,
    float phosphorus_affinity_p,
    void* params)
{
    Input* input = (Input*)params;
    float loc_phosphorus_nutrient_uptake = phosphorus_nutrient_uptake(input->max_phosphorus_uptake, input->phosphorus_soil_content, input->Michaelis_Menten_constant_phosphorus);
    float loc_pot_phosphorus_uptake = pot_phosphorus_uptake(loc_phosphorus_nutrient_uptake, phosphorus_affinity_p, input->root_biomass, input->total_biomass);
    float loc_phosphorus_uptake = phosphorus_uptake(loc_pot_phosphorus_uptake, input->sucrose, input->phosphorus_uptake_sucrose_consumption);

    float loc_min_phosphorus = min_phosphorus(
        input->respiration_frequency, input->max_sucrose, input->max_starch, input->sucrose_loading_frequency, input->assimilation_cost_phorphorus, input->min_phosphorus_photosynthesis, input->nutrient_conversion_parameter, input->photoperiod);

    float loc_max_phosphorus = max_phosphorus(input->max_nitrogen, input->optimal_stechiometric_ratio);

    return phosphorus_affinity(phosphorus_affinity_p, loc_phosphorus_uptake, input->phosphorus_cost, loc_max_phosphorus, input->phosphorus, input->photosynthesis, input->max_photosyntetic_rate, input->lambda_k, input->optimal_stechiometric_ratio, input->nitrogen, loc_min_phosphorus, input->phosphorus_uptake_sucrose_consumption, input->starch_partition_coeff, input->starch_degradation_rate);
}
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
    float min_nitrogen_photosyntesis)
{
    float first_term = (respiration_frequency * (sucrose + starch) + (night_efficiency_starch * sucrose_loading_frequency))
        * assimilation_cost_nitrogen * (leaf_biomass / (total_biomass + epsylon));

    float second_term = (photosynthesis / max_photosyntetic_rate) * min_nitrogen_photosyntesis;

    return first_term + second_term;
}

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
    float min_phosphorus_photosyntesis)
{
    float first_term = (respiration_frequency * (sucrose + starch) + (night_efficiency_starch * sucrose_loading_frequency))
        * assimilation_cost_phosphorus * (leaf_biomass / (total_biomass + epsylon));

    float second_term = (photosynthesis / max_photosyntetic_rate) * min_phosphorus_photosyntesis;

    return first_term + second_term;
}

float min_nitrogen(
    float respiration_frequency,
    float max_sucrose,
    float max_starch,
    float sucrose_loading_frequency,
    float assimilation_cost_nitrogen,
    float min_nitrogen_photosynthesis,
    float nutrient_conversion_parameter,
    float photoperiod)
{

    float n_min = ((respiration_frequency * (max_sucrose + max_starch)) + (sucrose_loading_frequency * max_sucrose)) * (((1.0f + 0.035f) * assimilation_cost_nitrogen) + (min_nitrogen_photosynthesis * nutrient_conversion_parameter * photoperiod));

    return n_min;
}

float min_phosphorus(
    float respiration_frequency,
    float max_sucrose,
    float max_starch,
    float sucrose_loading_frequency,
    float assimilation_cost_phosphorus,
    float min_phosphorus_photosynthesis,
    float nutrient_conversion_parameter,
    float photoperiod)
{

    float n_min = ((respiration_frequency * (max_sucrose + max_starch)) + (sucrose_loading_frequency * max_sucrose)) * (((1.0f + 0.035f) * assimilation_cost_phosphorus) + (min_phosphorus_photosynthesis * nutrient_conversion_parameter * photoperiod));

    return n_min;
}

float max_nitrogen(float min_nitrogen)
{
    return min_nitrogen * 4.0f;
}

float max_phosphorus(float max_nitrogen, float optimal_stechiometric_ratio)
{
    return max_nitrogen / optimal_stechiometric_ratio;
}

float starch_production(
    float photosynthesis,
    float starch_partition_coeff,
    float starch_degradation_rate)
{

    return (starch_partition_coeff * photosynthesis) - starch_degradation_rate;
}
float starch_production_f(
    float t,
    float starch,
    void* params)
{
    Input* input = (Input*)params;
    return starch_production(input->photosynthesis, input->starch_partition_coeff, input->starch_degradation_rate);
}

float sucrose_production(
    float sucrose_partition_coeff,
    float photosynthesis,
    float starch_degradation_rate,
    float uptake_cost,
    float transport_cost,
    float respiration_frequency,
    float sucrose,
    float sucrose_loading_frequency,
    float night_efficiency_starch)
{
    //printf(
        // "sucrose_partition_coeff=%8.3f  photosynthesis=%8.3f  starch_degradation_rate=%8.3f  "
        // "uptake_cost=%8.3f  transport_cost=%8.3f  respiration_frequency=%8.3f  "
        // "sucrose=%8.3f  sucrose_loading_frequency=%8.3f  night_efficiency_starch=%8.3f\n",
        // sucrose_partition_coeff,
        // photosynthesis,
        // starch_degradation_rate,
        // uptake_cost,
        // transport_cost,
        // respiration_frequency,
        // sucrose,
        // sucrose_loading_frequency,
        // night_efficiency_starch);

    return ((1.0f - sucrose_partition_coeff) * photosynthesis) + starch_degradation_rate - uptake_cost - transport_cost - (respiration_frequency * sucrose) - (sucrose_loading_frequency * night_efficiency_starch);
}
float sucrose_production_f(
    float t,
    float sucrose,
    void* params)
{
    Input* in = (Input*)params;    

    float loc_photosynthesis = photosynthesis(
        in->light,
        in->limitation_of_photosyntetic_rate,
        in->max_photosyntetic_rate,
        in->nitrogen_saturation,
        in->phosphorus_saturation,
        in->leaf_biomass,
        in->min_leaf_biomass);
    float loc_night_efficiency_starch = night_efficieny_starch(sucrose, in->max_sucrose, in->lambda_g, in->light, in->starch_partition_coeff);
    float loc_uptake_cost = uptake_cost(in->nitrogen_uptake_sucrose_consumption,
        in->nitrogen_uptake,
        in->phosphorus_uptake_sucrose_consumption,
        in->phosphorus_uptake);

    float loc_transport_cost = transport_cost(in->sucrose_consumption_transport,
        in->respiration_frequency,
        sucrose, /* <-- ovde! */
        in->sucrose_loading_frequency,
        loc_night_efficiency_starch);

    /* 2. Više NIŠTA ne upisujemo u in->...; derivacija mora biti čista */
    return sucrose_production(in->starch_partition_coeff,
        loc_photosynthesis,
        in->starch_degradation_rate,
        loc_uptake_cost,
        loc_transport_cost,
        in->respiration_frequency,
        sucrose,
        in->sucrose_loading_frequency,
        loc_night_efficiency_starch);
}

float night_efficieny_starch(
    float sucrose,
    float max_sucrose,
    float lambda_g,
    float light,
    float starch_partition_coeff)
{

    return fminf(sucrose, max_sucrose) * (lambda_g + ((1.0f - lambda_g) * (1.0f - light + (light * (1.0f - starch_partition_coeff)))));
}

float uptake_cost(
    float nitrogen_uptake_sucrose_consumption,
    float nitrogen_uptake,
    float phosphorus_uptake_sucrose_consumption,
    float phosphorus_uptake)
{
    return (nitrogen_uptake_sucrose_consumption * nitrogen_uptake) + (phosphorus_uptake_sucrose_consumption * phosphorus_uptake);
}

float transport_cost(
    float sucrose_consumption_transport,
    float respiration_frequency,
    float sucrose,
    float sucrose_loading_frequency,
    float night_efficiency_starch)
{
    return sucrose_consumption_transport * ((respiration_frequency * sucrose) + (sucrose_loading_frequency * night_efficiency_starch));
}

float leaf_growth(
    float labmda_sb,
    float sucrose_root_allocation,
    float sucrose_loading_frequency,
    float night_efficiency_starch,
    float leaf_biomass,
    float leaf_death_rate,
    float leaf_competative_rate)
{
        printf("\n[LEAF GROWTH DEBUG] lambda_sb: %.6f", labmda_sb);
    printf(" sucrose_root_allocation: %.6f", sucrose_root_allocation);
    printf(" sucrose_loading_frequency: %.6f", sucrose_loading_frequency);
    printf(" night_efficiency_starch: %.6f", night_efficiency_starch);
    printf(" leaf_biomass: %.6f", leaf_biomass);
    printf(" leaf_death_rate: %.6f", leaf_death_rate);
    printf(" leaf_competative_rate: %.6f", leaf_competative_rate);


    return (labmda_sb * (1.0f - sucrose_root_allocation) * sucrose_loading_frequency * night_efficiency_starch * leaf_biomass) - (leaf_death_rate * leaf_biomass) - (leaf_competative_rate * leaf_biomass * leaf_biomass);
}
float leaf_growth_f(
    float t,
    float leaf_biomass,
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

float root_growth(
    float labmda_sb,
    float sucrose_root_allocation,
    float sucrose_loading_frequency,
    float night_efficiency_starch,
    float leaf_biomass,
    float root_biomass,
    float root_death_rate,
    float root_competative_rate)
{

    printf(
        "\n[ROOT GROWTH DEBUG]"
        " lambda_sb= %.6f"
        " sucrose_root_allocation= %.6f"
        " sucrose_loading_frequency= %.6f"
        " night_efficiency_starch= %.6f"
        " leaf_biomass= %.6f"
        " root_biomass= %.6f"
        " root_death_rate= %.6f"
        " root_competative_rate= %.6f",
        labmda_sb,
        sucrose_root_allocation,
        sucrose_loading_frequency,
        night_efficiency_starch,
        leaf_biomass,
        root_biomass,
        root_death_rate,
        root_competative_rate
    );
    return (labmda_sb * sucrose_root_allocation * sucrose_loading_frequency * night_efficiency_starch * leaf_biomass) - (root_death_rate * root_biomass) - (root_competative_rate * root_biomass * root_biomass);
}
float root_growth_f(
    float t,
    float root_biomass,
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

float stochiometric_signal(
    float optimal_stechiometric_ratio,
    float nitrogen,
    float phosphorus)
{
    return optimal_stechiometric_ratio / (optimal_stechiometric_ratio + (nitrogen / (phosphorus + epsylon)));
}

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
    float min_phosphorus)
{
    printf(
        "\n[SUCROSE ROOT ALLOCATION] sucrose_root_allocation=%8.3f  nitrogen_affinity=%8.3f  phospsorus_affinity=%8.3f  "
        "stochiometric_signal=%8.3f  sucrose=%8.3f  min_sucrose=%8.3f  "
        "nitrogen=%8.3f  min_nitrogen=%8.3f  phosphorus=%8.3f  min_phosphorus=%8.3f\n",
        sucrose_root_allocation,
        nitrogen_affinity,
        phospsorus_affinity,
        stochiometric_signal,
        sucrose,
        min_sucrose,
        nitrogen,
        min_nitrogen,
        phosphorus,
        min_phosphorus);

    float first = (1.0f - sucrose_root_allocation) * ((nitrogen_affinity * stochiometric_signal) + ((1.0f - stochiometric_signal) * phospsorus_affinity));
    float secound = sucrose_root_allocation * (((nitrogen * stochiometric_signal) / (nitrogen + min_nitrogen)) + ((phosphorus * (1.0f - stochiometric_signal)) / (phosphorus + min_phosphorus)) + (min_sucrose / (min_sucrose + sucrose)));
    return first - secound;
}
float sucrose_root_allocation_f(
    float t,
    float sucrose_root_allocation_p,
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

// Farquhar C3 Photosynthesis Model


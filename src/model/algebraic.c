#include "model/algebraic.h"

#include "model/model.h"

static Input plant_input_from_typed_data(
    const PlantState* state,
    const PlantParameters* parameters)
{
    Input input = {
        .core = parameters->core,
        .gas_exchange = parameters->gas_exchange
    };
    plant_state_apply_to_input(&input, state);
    return input;
}

void plant_compute_outputs(
    const PlantState* state,
    const PlantParameters* parameters,
    const PlantEnvironment* environment,
    real_t starch_night_start,
    PlantOutputs* outputs)
{
    Input s = plant_input_from_typed_data(state, parameters);

    update_light_conditions(&s, environment->hour_of_day);

    s.core.min_nitrogen = min_nitrogen(
        s.core.respiration_frequency, s.core.max_sucrose, s.core.max_starch,
        s.core.sucrose_loading_frequency, s.core.assimilation_cost_nitrogen,
        s.core.min_nitrogen_photosynthesis, s.core.nutrient_conversion_parameter, s.core.photoperiod);

    s.core.min_phosphorus = min_phosphorus(
        s.core.respiration_frequency, s.core.max_sucrose, s.core.max_starch,
        s.core.sucrose_loading_frequency, s.core.assimilation_cost_phosphorus,
        s.core.min_phosphorus_photosynthesis, s.core.nutrient_conversion_parameter, s.core.photoperiod);

    s.core.max_starch = max_starch(s.core.max_starch_degradation_rate, s.core.photoperiod);
    s.core.max_nitrogen = max_nitrogen(s.core.min_nitrogen);
    s.core.max_phosphorus = max_phosphorus(s.core.max_nitrogen, s.core.optimal_stoichiometric_ratio);

    s.core.nitrogen_saturation = nitrogen_saturation(s.core.nitrogen, s.core.min_nitrogen_photosynthesis);
    s.core.phosphorus_saturation = phosphorus_saturation(
        s.core.min_nitrogen_photosynthesis, s.core.optimal_stoichiometric_ratio, s.core.phosphorus);

    s.core.limitation_of_photosynthetic_rate = limitation_of_photosynthetic_rate(s.core.starch, s.core.feedback_on_photosynthesis, s.core.max_starch);

    s.core.photosynthesis = farquhar_photosynthesis(&s);

    s.core.night_efficiency_starch = night_efficiency_starch(
        s.core.sucrose, s.core.max_sucrose, s.core.lambda_g, s.core.light, s.core.starch_partition_coeff);

    s.core.starch_degradation_rate = starch_degradation(
        s.core.max_starch_degradation_rate,
        s.core.max_sucrose,
        starch_night_start,
        s.core.min_starch,
        s.core.sucrose,
        s.core.light,
        environment->hour_of_day,
        s.core.photoperiod);

    s.core.uptake_cost = uptake_cost(
        s.core.nitrogen_uptake_sucrose_consumption, s.core.nitrogen_uptake,
        s.core.phosphorus_uptake_sucrose_consumption, s.core.phosphorus_uptake);

    s.core.transport_cost = transport_cost(
        s.core.sucrose_consumption_transport, s.core.respiration_frequency,
        s.core.sucrose, s.core.sucrose_loading_frequency, s.core.night_efficiency_starch);

    s.core.nitrogen_cost = nitrogen_cost(
        s.core.respiration_frequency, s.core.sucrose, s.core.starch,
        s.core.sucrose_loading_frequency, s.core.night_efficiency_starch,
        s.core.assimilation_cost_nitrogen,
        s.core.leaf_biomass, s.core.total_biomass,
        s.core.photosynthesis, s.core.max_photosynthetic_rate,
        s.core.min_nitrogen_photosynthesis);

    s.core.phosphorus_cost = phosphorus_cost(
        s.core.respiration_frequency, s.core.sucrose, s.core.starch,
        s.core.sucrose_loading_frequency, s.core.night_efficiency_starch,
        s.core.assimilation_cost_phosphorus,
        s.core.leaf_biomass, s.core.total_biomass,
        s.core.photosynthesis, s.core.max_photosynthetic_rate,
        s.core.min_phosphorus_photosynthesis);

    s.core.nitrogen_nutrient_uptake = nitrogen_nutrient_uptake(
        s.core.max_nitrogen_uptake, s.core.nitrogen_soil_content, s.core.Michaelis_Menten_constant_nitrogen);

    s.core.phosphorus_nutrient_uptake = phosphorus_nutrient_uptake(
        s.core.max_phosphorus_uptake, s.core.phosphorus_soil_content, s.core.Michaelis_Menten_constant_phosphorus);

    s.core.pot_nitrogen_uptake = pot_nitrogen_uptake(
        s.core.nitrogen_nutrient_uptake, s.core.nitrogen_affinity, s.core.root_biomass, s.core.total_biomass);

    s.core.pot_phosphorus_uptake = pot_phosphorus_uptake(
        s.core.phosphorus_nutrient_uptake, s.core.phosphorus_affinity, s.core.root_biomass, s.core.total_biomass);

    s.core.nitrogen_uptake = nitrogen_uptake(
        s.core.pot_nitrogen_uptake, s.core.sucrose, s.core.nitrogen_uptake_sucrose_consumption);

    s.core.phosphorus_uptake = phosphorus_uptake(
        s.core.pot_phosphorus_uptake, s.core.sucrose, s.core.phosphorus_uptake_sucrose_consumption);

    s.core.stoichiometric_signal = stoichiometric_signal(
        s.core.optimal_stoichiometric_ratio, s.core.nitrogen, s.core.phosphorus);

    outputs->light = s.core.light;
    outputs->light_PAR = s.core.light_PAR;
    outputs->limitation_of_photosynthetic_rate = s.core.limitation_of_photosynthetic_rate;
    outputs->max_starch = s.core.max_starch;
    outputs->max_nitrogen = s.core.max_nitrogen;
    outputs->min_nitrogen = s.core.min_nitrogen;
    outputs->max_phosphorus = s.core.max_phosphorus;
    outputs->min_phosphorus = s.core.min_phosphorus;
    outputs->nitrogen_saturation = s.core.nitrogen_saturation;
    outputs->phosphorus_saturation = s.core.phosphorus_saturation;
    outputs->photosynthesis = s.core.photosynthesis;
    outputs->night_efficiency_starch = s.core.night_efficiency_starch;
    outputs->starch_degradation_rate = s.core.starch_degradation_rate;
    outputs->uptake_cost = s.core.uptake_cost;
    outputs->transport_cost = s.core.transport_cost;
    outputs->nitrogen_cost = s.core.nitrogen_cost;
    outputs->phosphorus_cost = s.core.phosphorus_cost;
    outputs->nitrogen_nutrient_uptake = s.core.nitrogen_nutrient_uptake;
    outputs->phosphorus_nutrient_uptake = s.core.phosphorus_nutrient_uptake;
    outputs->pot_nitrogen_uptake = s.core.pot_nitrogen_uptake;
    outputs->pot_phosphorus_uptake = s.core.pot_phosphorus_uptake;
    outputs->nitrogen_uptake = s.core.nitrogen_uptake;
    outputs->phosphorus_uptake = s.core.phosphorus_uptake;
    outputs->stoichiometric_signal = s.core.stoichiometric_signal;
    outputs->total_biomass = s.core.total_biomass;
}

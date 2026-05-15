#include "model/algebraic.h"

#include "model/model.h"

static Input plant_input_from_typed_data(
    const PlantState* state,
    const PlantParameters* parameters)
{
    Input input = {
        .photo = parameters->photo,
        .carbohydrates = parameters->carbohydrates,
        .nutrients = parameters->nutrients,
        .growth = parameters->growth,
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

    s.nutrients.min_nitrogen = min_nitrogen(
        s.carbohydrates.respiration_frequency, s.carbohydrates.max_sucrose, s.carbohydrates.max_starch,
        s.carbohydrates.sucrose_loading_frequency, s.nutrients.assimilation_cost_nitrogen,
        s.nutrients.min_nitrogen_photosynthesis, s.nutrients.nutrient_conversion_parameter, s.photo.photoperiod);

    s.nutrients.min_phosphorus = min_phosphorus(
        s.carbohydrates.respiration_frequency, s.carbohydrates.max_sucrose, s.carbohydrates.max_starch,
        s.carbohydrates.sucrose_loading_frequency, s.nutrients.assimilation_cost_phosphorus,
        s.nutrients.min_phosphorus_photosynthesis, s.nutrients.nutrient_conversion_parameter, s.photo.photoperiod);

    s.carbohydrates.max_starch = max_starch(s.carbohydrates.max_starch_degradation_rate, s.photo.photoperiod);
    s.nutrients.max_nitrogen = max_nitrogen(s.nutrients.min_nitrogen);
    s.nutrients.max_phosphorus = max_phosphorus(s.nutrients.max_nitrogen, s.nutrients.optimal_stoichiometric_ratio);

    s.photo.nitrogen_saturation = nitrogen_saturation(s.nutrients.nitrogen, s.nutrients.min_nitrogen_photosynthesis);
    s.photo.phosphorus_saturation = phosphorus_saturation(
        s.nutrients.min_nitrogen_photosynthesis, s.nutrients.optimal_stoichiometric_ratio, s.nutrients.phosphorus);

    s.photo.limitation_of_photosynthetic_rate = limitation_of_photosynthetic_rate(s.carbohydrates.starch, s.photo.feedback_on_photosynthesis, s.carbohydrates.max_starch);

    s.photo.photosynthesis = farquhar_photosynthesis(&s);

    s.carbohydrates.night_efficiency_starch = night_efficiency_starch(
        s.carbohydrates.sucrose, s.carbohydrates.max_sucrose, s.carbohydrates.lambda_g, s.photo.light, s.carbohydrates.starch_partition_coeff);

    s.carbohydrates.starch_degradation_rate = starch_degradation(
        s.carbohydrates.max_starch_degradation_rate,
        s.carbohydrates.max_sucrose,
        starch_night_start,
        s.carbohydrates.min_starch,
        s.carbohydrates.sucrose,
        s.photo.light,
        environment->hour_of_day,
        s.photo.photoperiod);

    s.carbohydrates.uptake_cost = uptake_cost(
        s.nutrients.nitrogen_uptake_sucrose_consumption, s.nutrients.nitrogen_uptake,
        s.nutrients.phosphorus_uptake_sucrose_consumption, s.nutrients.phosphorus_uptake);

    s.carbohydrates.transport_cost = transport_cost(
        s.carbohydrates.sucrose_consumption_transport, s.carbohydrates.respiration_frequency,
        s.carbohydrates.sucrose, s.carbohydrates.sucrose_loading_frequency, s.carbohydrates.night_efficiency_starch);

    s.nutrients.nitrogen_cost = nitrogen_cost(
        s.carbohydrates.respiration_frequency, s.carbohydrates.sucrose, s.carbohydrates.starch,
        s.carbohydrates.sucrose_loading_frequency, s.carbohydrates.night_efficiency_starch,
        s.nutrients.assimilation_cost_nitrogen,
        s.growth.leaf_biomass, s.growth.total_biomass,
        s.photo.photosynthesis, s.photo.max_photosynthetic_rate,
        s.nutrients.min_nitrogen_photosynthesis);

    s.nutrients.phosphorus_cost = phosphorus_cost(
        s.carbohydrates.respiration_frequency, s.carbohydrates.sucrose, s.carbohydrates.starch,
        s.carbohydrates.sucrose_loading_frequency, s.carbohydrates.night_efficiency_starch,
        s.nutrients.assimilation_cost_phosphorus,
        s.growth.leaf_biomass, s.growth.total_biomass,
        s.photo.photosynthesis, s.photo.max_photosynthetic_rate,
        s.nutrients.min_phosphorus_photosynthesis);

    s.nutrients.nitrogen_nutrient_uptake = nitrogen_nutrient_uptake(
        s.nutrients.max_nitrogen_uptake, s.nutrients.nitrogen_soil_content, s.nutrients.Michaelis_Menten_constant_nitrogen);

    s.nutrients.phosphorus_nutrient_uptake = phosphorus_nutrient_uptake(
        s.nutrients.max_phosphorus_uptake, s.nutrients.phosphorus_soil_content, s.nutrients.Michaelis_Menten_constant_phosphorus);

    s.nutrients.pot_nitrogen_uptake = pot_nitrogen_uptake(
        s.nutrients.nitrogen_nutrient_uptake, s.nutrients.nitrogen_affinity, s.growth.root_biomass, s.growth.total_biomass);

    s.nutrients.pot_phosphorus_uptake = pot_phosphorus_uptake(
        s.nutrients.phosphorus_nutrient_uptake, s.nutrients.phosphorus_affinity, s.growth.root_biomass, s.growth.total_biomass);

    s.nutrients.nitrogen_uptake = nitrogen_uptake(
        s.nutrients.pot_nitrogen_uptake, s.carbohydrates.sucrose, s.nutrients.nitrogen_uptake_sucrose_consumption);

    s.nutrients.phosphorus_uptake = phosphorus_uptake(
        s.nutrients.pot_phosphorus_uptake, s.carbohydrates.sucrose, s.nutrients.phosphorus_uptake_sucrose_consumption);

    s.nutrients.stoichiometric_signal = stoichiometric_signal(
        s.nutrients.optimal_stoichiometric_ratio, s.nutrients.nitrogen, s.nutrients.phosphorus);

    outputs->light = s.photo.light;
    outputs->light_PAR = s.photo.light_PAR;
    outputs->limitation_of_photosynthetic_rate = s.photo.limitation_of_photosynthetic_rate;
    outputs->max_starch = s.carbohydrates.max_starch;
    outputs->max_nitrogen = s.nutrients.max_nitrogen;
    outputs->min_nitrogen = s.nutrients.min_nitrogen;
    outputs->max_phosphorus = s.nutrients.max_phosphorus;
    outputs->min_phosphorus = s.nutrients.min_phosphorus;
    outputs->nitrogen_saturation = s.photo.nitrogen_saturation;
    outputs->phosphorus_saturation = s.photo.phosphorus_saturation;
    outputs->photosynthesis = s.photo.photosynthesis;
    outputs->night_efficiency_starch = s.carbohydrates.night_efficiency_starch;
    outputs->starch_degradation_rate = s.carbohydrates.starch_degradation_rate;
    outputs->uptake_cost = s.carbohydrates.uptake_cost;
    outputs->transport_cost = s.carbohydrates.transport_cost;
    outputs->nitrogen_cost = s.nutrients.nitrogen_cost;
    outputs->phosphorus_cost = s.nutrients.phosphorus_cost;
    outputs->nitrogen_nutrient_uptake = s.nutrients.nitrogen_nutrient_uptake;
    outputs->phosphorus_nutrient_uptake = s.nutrients.phosphorus_nutrient_uptake;
    outputs->pot_nitrogen_uptake = s.nutrients.pot_nitrogen_uptake;
    outputs->pot_phosphorus_uptake = s.nutrients.pot_phosphorus_uptake;
    outputs->nitrogen_uptake = s.nutrients.nitrogen_uptake;
    outputs->phosphorus_uptake = s.nutrients.phosphorus_uptake;
    outputs->stoichiometric_signal = s.nutrients.stoichiometric_signal;
    outputs->total_biomass = s.growth.total_biomass;
}

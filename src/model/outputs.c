#include "model/outputs.h"

void plant_outputs_apply_to_input(Input* input, const PlantOutputs* outputs)
{
    input->photo.light = outputs->light;
    input->photo.light_PAR = outputs->light_PAR;
    input->photo.limitation_of_photosynthetic_rate = outputs->limitation_of_photosynthetic_rate;
    input->carbohydrates.max_starch = outputs->max_starch;
    input->nutrients.max_nitrogen = outputs->max_nitrogen;
    input->nutrients.min_nitrogen = outputs->min_nitrogen;
    input->nutrients.max_phosphorus = outputs->max_phosphorus;
    input->nutrients.min_phosphorus = outputs->min_phosphorus;
    input->photo.nitrogen_saturation = outputs->nitrogen_saturation;
    input->photo.phosphorus_saturation = outputs->phosphorus_saturation;
    input->photo.photosynthesis = outputs->photosynthesis;
    input->carbohydrates.night_efficiency_starch = outputs->night_efficiency_starch;
    input->carbohydrates.starch_degradation_rate = outputs->starch_degradation_rate;
    input->carbohydrates.uptake_cost = outputs->uptake_cost;
    input->carbohydrates.transport_cost = outputs->transport_cost;
    input->nutrients.nitrogen_cost = outputs->nitrogen_cost;
    input->nutrients.phosphorus_cost = outputs->phosphorus_cost;
    input->nutrients.nitrogen_nutrient_uptake = outputs->nitrogen_nutrient_uptake;
    input->nutrients.phosphorus_nutrient_uptake = outputs->phosphorus_nutrient_uptake;
    input->nutrients.pot_nitrogen_uptake = outputs->pot_nitrogen_uptake;
    input->nutrients.pot_phosphorus_uptake = outputs->pot_phosphorus_uptake;
    input->nutrients.nitrogen_uptake = outputs->nitrogen_uptake;
    input->nutrients.phosphorus_uptake = outputs->phosphorus_uptake;
    input->nutrients.stoichiometric_signal = outputs->stoichiometric_signal;
    input->growth.total_biomass = outputs->total_biomass;
}

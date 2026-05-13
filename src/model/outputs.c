#include "model/outputs.h"

void plant_outputs_apply_to_input(Input* input, const PlantOutputs* outputs)
{
    input->core.light = outputs->light;
    input->core.light_PAR = outputs->light_PAR;
    input->core.limitation_of_photosynthetic_rate = outputs->limitation_of_photosynthetic_rate;
    input->core.max_starch = outputs->max_starch;
    input->core.max_nitrogen = outputs->max_nitrogen;
    input->core.min_nitrogen = outputs->min_nitrogen;
    input->core.max_phosphorus = outputs->max_phosphorus;
    input->core.min_phosphorus = outputs->min_phosphorus;
    input->core.nitrogen_saturation = outputs->nitrogen_saturation;
    input->core.phosphorus_saturation = outputs->phosphorus_saturation;
    input->core.photosynthesis = outputs->photosynthesis;
    input->core.night_efficiency_starch = outputs->night_efficiency_starch;
    input->core.starch_degradation_rate = outputs->starch_degradation_rate;
    input->core.uptake_cost = outputs->uptake_cost;
    input->core.transport_cost = outputs->transport_cost;
    input->core.nitrogen_cost = outputs->nitrogen_cost;
    input->core.phosphorus_cost = outputs->phosphorus_cost;
    input->core.nitrogen_nutrient_uptake = outputs->nitrogen_nutrient_uptake;
    input->core.phosphorus_nutrient_uptake = outputs->phosphorus_nutrient_uptake;
    input->core.pot_nitrogen_uptake = outputs->pot_nitrogen_uptake;
    input->core.pot_phosphorus_uptake = outputs->pot_phosphorus_uptake;
    input->core.nitrogen_uptake = outputs->nitrogen_uptake;
    input->core.phosphorus_uptake = outputs->phosphorus_uptake;
    input->core.stoichiometric_signal = outputs->stoichiometric_signal;
    input->core.total_biomass = outputs->total_biomass;
}

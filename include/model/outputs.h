#ifndef MODEL_OUTPUTS_H
#define MODEL_OUTPUTS_H

#include "model/input.h"

typedef struct PlantOutputs {
    real_t light;
    real_t light_PAR;
    real_t limitation_of_photosynthetic_rate;
    real_t max_starch;
    real_t max_nitrogen;
    real_t min_nitrogen;
    real_t max_phosphorus;
    real_t min_phosphorus;
    real_t nitrogen_saturation;
    real_t phosphorus_saturation;
    real_t photosynthesis;
    real_t night_efficiency_starch;
    real_t starch_degradation_rate;
    real_t uptake_cost;
    real_t transport_cost;
    real_t nitrogen_cost;
    real_t phosphorus_cost;
    real_t nitrogen_nutrient_uptake;
    real_t phosphorus_nutrient_uptake;
    real_t pot_nitrogen_uptake;
    real_t pot_phosphorus_uptake;
    real_t nitrogen_uptake;
    real_t phosphorus_uptake;
    real_t stoichiometric_signal;
    real_t total_biomass;
} PlantOutputs;

void plant_outputs_apply_to_input(Input* input, const PlantOutputs* outputs);

#endif

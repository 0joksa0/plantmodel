#ifndef MODEL_ALGEBRAIC_H
#define MODEL_ALGEBRAIC_H

#include "model/environment.h"
#include "model/outputs.h"
#include "model/parameters.h"
#include "model/state.h"

void plant_compute_outputs(
    const PlantState* state,
    const PlantParameters* parameters,
    const PlantEnvironment* environment,
    real_t starch_night_start,
    PlantOutputs* outputs);

#endif

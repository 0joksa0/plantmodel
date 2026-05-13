#ifndef MODEL_ENVIRONMENT_H
#define MODEL_ENVIRONMENT_H

#include "model/input.h"
#include <solver.h>

typedef struct PlantEnvironment {
    real_t hour_of_day;
} PlantEnvironment;

PlantEnvironment plant_environment_from_input(const Input* input, real_t hour_of_day);

#endif

#ifndef MODEL_PARAMETERS_H
#define MODEL_PARAMETERS_H

#include "model/input.h"

typedef struct PlantParameters {
    InputCore core;
    InputGasExchange gas_exchange;
} PlantParameters;

PlantParameters plant_parameters_from_input(const Input* input);

#endif

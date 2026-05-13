#include "model/parameters.h"

PlantParameters plant_parameters_from_input(const Input* input)
{
    PlantParameters parameters = {
        .core = input->core,
        .gas_exchange = input->gas_exchange
    };
    return parameters;
}

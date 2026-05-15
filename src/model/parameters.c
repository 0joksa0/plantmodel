#include "model/parameters.h"

PlantParameters plant_parameters_from_input(const Input* input)
{
    PlantParameters parameters = {
        .photo = input->photo,
        .carbohydrates = input->carbohydrates,
        .nutrients = input->nutrients,
        .growth = input->growth,
        .gas_exchange = input->gas_exchange
    };
    return parameters;
}

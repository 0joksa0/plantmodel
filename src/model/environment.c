#include "model/environment.h"

PlantEnvironment plant_environment_from_input(const Input* input, real_t hour_of_day)
{
    (void)input;
    PlantEnvironment environment = {
        .hour_of_day = hour_of_day
    };
    return environment;
}

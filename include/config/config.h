#ifndef PROJECT_CONFIG_H
#define PROJECT_CONFIG_H

#include "model/input.h"
#include "simulation/simulation.h"
#include <stddef.h>

int config_load_file(const char* path, SimulationConfig* sim, Input* input);
int config_get_input_field_offset(const char* field_name, size_t* offset_out);

#endif

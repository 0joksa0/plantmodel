#define SOLVER_REAL_AS double
#include "config/config.h"
#include "model/input.h"
#include "model/model.h"
#include "simulation/simulation.h"
#include <solver.h>
#include <stdio.h>
#include <string.h>

int main(int argc, char** argv)
{
    SimulationConfig config = simulation_default_config();
    Input input = generate_input();
    const char* config_path = "config/default.config";

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--config") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Missing value for --config\n");
                return 1;
            }
            config_path = argv[++i];
        }
    }

    if (config_load_file(config_path, &config, &input) != 0) {
        fprintf(stderr, "Failed to load config: %s\n", config_path);
        return 1;
    }

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--no-gui") == 0) {
            config.gui_enabled = 0;
        } else if (strcmp(argv[i], "--gui") == 0) {
            config.gui_enabled = 1;
        } else if (strcmp(argv[i], "--config") == 0) {
            ++i;
        } else if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
            printf("Usage: %s [--config <path>] [--gui | --no-gui]\n", argv[0]);
            return 0;
        }
    }

    printf("=== Simulacija rasta biljke ===\n");
    photo_dataset_open("photo_surrogate_dataset.csv");

    real_t FW_start = input.core.total_biomass;

    SimulationResult result;
    simulate_days(&config, &input, &result);

    real_t duration_days = (real_t)config.days;
    real_t rgr = compute_rgr(FW_start, input.core.total_biomass, duration_days);
    printf("%.1f model_RGR: %.6f, leaf_start=%.6f, leaf_end=%.6f\n",
        (double)input.core.photoperiod, (double)rgr, (double)FW_start, (double)input.core.total_biomass);

    photo_dataset_close();
    simulation_result_free(&result);
    return 0;
}

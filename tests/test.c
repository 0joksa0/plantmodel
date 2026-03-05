#include <criterion/criterion.h>
#include <solver.h>
#include <stdio.h>

#include "model/input.h"
#include "simulation/simulation.h"

#define TEST_RGR(hours, lambda_sni_v, lambda_sb_v, feedback_v, lit_val)                \
    Test(math, hours##_hours)                                                          \
    {                                                                                  \
        /* TODO: Add assertions for final state ranges (N/P/sucrose/starch) to catch regressions not visible in RGR alone. */ \
        printf(#hours "h\n");                                                          \
        Input input = generate_input();                                                \
        input.core.photoperiod = REAL(hours);                                               \
        input.core.lambda_sni = REAL(lambda_sni_v);                                         \
        input.core.lambda_sb = REAL(lambda_sb_v);                                           \
        input.core.feedback_on_photosynthesis = REAL(feedback_v);                           \
        SimulationConfig config = simulation_default_config();                          \
        config.days = 30;                                                               \
        config.gui_enabled = 0;                                                         \
        real_t FW_start = input.core.total_biomass;                                          \
        SimulationResult result;                                                        \
        simulate_days(&config, &input, &result);                                        \
        real_t rgr = compute_rgr(FW_start, input.core.total_biomass, 29.0);                 \
        printf("%d model_RGR: %.6f, lit_RGR: %.6f\n",hours, rgr, REAL(lit_val));                \
        printf("%d leaf_start=%.6f, leaf_end=%.6f\n",hours, FW_start, input.core.total_biomass); \
        simulation_result_free(&result);                                               \
        cr_assert(fabs(rgr - REAL(lit_val)) < lit_val / REAL(10),                                    \
            "RGR mismatch at " #hours "h: model=%.6f, lit=%.6f", rgr, REAL(lit_val));  \
    }

TEST_RGR(4, 0.16, 0.00344, 0.82, 0.0680)
TEST_RGR(6, 0.15, 0.00413, 0.79, 0.1135)
TEST_RGR(8, 0.13, 0.00515, 0.67, 0.1708)
TEST_RGR(12, 0.08, 0.00578, 0.62, 0.2600)
// TEST_RGR(18, 0.004, 0.00524, 0.55, 0.3065)

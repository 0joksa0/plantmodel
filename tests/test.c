#include <criterion/criterion.h>
#include <solver.h>
#include <stdio.h>

#include "model/input.h"

#define PHOTO_PERIOD REAL(8.0)
void simulate_days(int days, Input* input);
real_t compute_rgr(real_t FW_start, real_t FW_end, real_t duration_hours)
{
    return (RLOG(FW_end) - RLOG(FW_start)) / duration_hours;
}

#define TEST_RGR(hours, lambda_sni_v, lambda_sb_v, feedback_v, lit_val)                \
    Test(math, hours##_hours)                                                          \
    {                                                                                  \
        printf(#hours "h\n");                                                          \
        Input input = generate_input();                                                \
        input.photoperiod = REAL(hours);                                               \
        input.lambda_sni = REAL(lambda_sni_v);                                         \
        input.lambda_sb = REAL(lambda_sb_v);                                           \
        input.feedback_on_photosynthesis = REAL(feedback_v);                           \
        real_t FW_start = input.total_biomass;                                          \
        simulate_days(30, &input);                                                     \
        real_t rgr = compute_rgr(FW_start, input.total_biomass, 30.0);                 \
        printf("model_RGR: %.6f, lit_RGR: %.6f\n", rgr, REAL(lit_val));                \
        printf("leaf_start=%.6f, leaf_end=%.6f\n", FW_start, input.total_biomass); \
        cr_assert(fabs(rgr - REAL(lit_val)) < 0.01,                                    \
            "RGR mismatch at " #hours "h: model=%.6f, lit=%.6f", rgr, REAL(lit_val));  \
    }

TEST_RGR(4, 0.16, 0.00344, 0.82, 0.0680)
TEST_RGR(6, 0.15, 0.00413, 0.79, 0.1135)
TEST_RGR(8, 0.13, 0.00515, 0.67, 0.1708)
TEST_RGR(12, 0.08, 0.00578, 0.62, 0.2600)
TEST_RGR(18, 0.004, 0.00524, 0.55, 0.3065)

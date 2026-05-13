#include "model/photo_surrogate.h"
#include "model/photo_surrogate_export.h"

#include <math.h>

#define PHOTO_SURROGATE_MIN_LIGHT_PAR 70.0

typedef char photo_surrogate_requires_3_inputs[
    (PHOTO_HIGH_MODEL_N_INPUTS == 3) ? 1 : -1];
typedef char photo_surrogate_requires_1_output[
    (PHOTO_HIGH_MODEL_N_OUTPUTS == 1) ? 1 : -1];

double photo_surrogate_should_use(double light_PAR)
{
    return light_PAR >= PHOTO_SURROGATE_MIN_LIGHT_PAR ? 1.0 : 0.0;
}

double photo_surrogate_predict_an(
    double light_PAR,
    double leaf_biomass,
    double stomatal_conductance)
{
    double x[PHOTO_HIGH_MODEL_N_INPUTS];
    double hidden[PHOTO_HIGH_MODEL_N_H1];
    double y_scaled;
    int i;
    int j;

    x[0] = (light_PAR - photo_high_model_x_mean[0]) /
           photo_high_model_x_scale[0];
    x[1] = (leaf_biomass - photo_high_model_x_mean[1]) /
           photo_high_model_x_scale[1];
    x[2] = (stomatal_conductance - photo_high_model_x_mean[2]) /
           photo_high_model_x_scale[2];

    for (j = 0; j < PHOTO_HIGH_MODEL_N_H1; ++j) {
        double z = photo_high_model_b1[j];

        for (i = 0; i < PHOTO_HIGH_MODEL_N_INPUTS; ++i) {
            z += x[i] * photo_high_model_W1[i][j];
        }

        hidden[j] = tanh(z);
    }

    y_scaled = photo_high_model_b2[0];
    for (j = 0; j < PHOTO_HIGH_MODEL_N_H1; ++j) {
        y_scaled += hidden[j] * photo_high_model_W2[j][0];
    }

    return (y_scaled * photo_high_model_y_scale[0]) +
           photo_high_model_y_mean[0];
}

double photo_surrogate_predict_if_applicable(
    double light_PAR,
    double leaf_biomass,
    double stomatal_conductance,
    int *used_surrogate)
{
    if (!photo_surrogate_should_use(light_PAR)) {
        if (used_surrogate) {
            *used_surrogate = 0;
        }
        return 0.0;
    }

    if (used_surrogate) {
        *used_surrogate = 1;
    }

    return photo_surrogate_predict_an(
        light_PAR,
        leaf_biomass,
        stomatal_conductance);
}

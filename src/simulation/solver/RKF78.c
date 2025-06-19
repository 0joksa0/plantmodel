
#include "simulation/solver/RKF78.h"
#include <stdio.h>

#define STAGES 13

static const float c[STAGES] = {
    0.0f,
    2.0f / 27.0f,
    1.0f / 9.0f,
    1.0f / 6.0f,
    5.0f / 12.0f,
    0.5f,
    5.0f / 6.0f,
    1.0f,
    1.0f,
    1.0f,
    1.0f,
    1.0f,
    1.0f
};

static const float a[STAGES][STAGES] = {
    { 0 },
    { 2.0f / 27.0f },
    { 1.0f / 36.0f, 1.0f / 12.0f },
    { 1.0f / 24.0f, 0.0f, 1.0f / 8.0f },
    { 5.0f / 12.0f, 0.0f, -25.0f / 16.0f, 25.0f / 16.0f },
    { 1.0f / 20.0f, 0.0f, 0.0f, 1.0f / 4.0f, 1.0f / 5.0f },
    { -25.0f / 108.0f, 0.0f, 0.0f, 125.0f / 108.0f, -65.0f / 27.0f, 125.0f / 54.0f },
    { 31.0f / 300.0f, 0.0f, 0.0f, 0.0f, 61.0f / 225.0f, -2.0f / 9.0f, 13.0f / 900.0f },
    { 2.0f, 0.0f, 0.0f, -53.0f / 6.0f, 704.0f / 45.0f, -107.0f / 9.0f, 67.0f / 90.0f, 3.0f },
    { -91.0f / 108.0f, 0.0f, 0.0f, 23.0f / 108.0f, -976.0f / 135.0f, 311.0f / 54.0f, -19.0f / 60.0f, 17.0f / 6.0f, -1.0f / 12.0f },
    { 2383.0f / 4100.0f, 0.0f, 0.0f, -341.0f / 164.0f, 4496.0f / 1025.0f, -301.0f / 82.0f, 2133.0f / 4100.0f, 45.0f / 82.0f, 45.0f / 164.0f, 18.0f / 41.0f },
    { 3.0f / 205.0f, 0.0f, 0.0f, 0.0f, 0.0f, -6.0f / 41.0f, -3.0f / 205.0f, -3.0f / 41.0f, 3.0f / 41.0f, 6.0f / 41.0f, 0.0f },
    { -1777.0f / 4100.0f, 0.0f, 0.0f, -341.0f / 164.0f, 4496.0f / 1025.0f, -289.0f / 82.0f, 2193.0f / 4100.0f, 51.0f / 82.0f, 33.0f / 164.0f, 12.0f / 41.0f, 0.0f, 1.0f }
};

static const float b7[STAGES] = {
    41.0f / 840.0f, 0.0f, 0.0f, 0.0f, 0.0f, 34.0f / 105.0f, 9.0f / 35.0f,
    9.0f / 35.0f, 9.0f / 280.0f, 9.0f / 280.0f, 41.0f / 840.0f, 0.0f, 0.0f
};

static const float b8[STAGES] = {
    0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 34.0f / 105.0f, 9.0f / 35.0f,
    9.0f / 35.0f, 9.0f / 280.0f, 9.0f / 280.0f, 0.0f, 0.0f, 0.0f
};

float rkf78_solver(
    float y0, float t0, float t_end,
    float h_init, float tol,
    ODEFunction f, void* params)
{
    float t = t0;
    float y = y0;
    float h = h_init;

    const float safety = 0.9f;
    const float min_scale = 0.1f;
    const float max_scale = 4.0f;

    while (t < t_end) {

        if (t + h > t_end)
            h = t_end - t;

        float k[STAGES];
        k[0] = h * f(t, y, params);

        for (int i = 1; i < STAGES; ++i) {
            float sum = 0.0f;
            for (int j = 0; j < i; ++j) {
                sum += a[i][j] * k[j];
            }
            k[i] = h * f(t + c[i] * h, y + sum, params);
        }

        float y7 = y;
        float y8 = y;
        for (int i = 0; i < STAGES; ++i) {
            y7 += b7[i] * k[i];
            y8 += b8[i] * k[i];
        }

        float err = fabsf(y8 - y7);


        if (err < tol || h < 1e-6f) {
            t += h;
            y = y8;
        } else {
        }

        float scale = safety * powf(tol / (err + FLT_EPSILON), 1.0f / 8.0f);
        scale = fmaxf(min_scale, fminf(scale, max_scale));
        h *= scale;
    }

    return y;
}

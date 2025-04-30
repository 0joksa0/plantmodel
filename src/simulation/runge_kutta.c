#include "simulation/runge_cutta.h"

float runge_kutta_4(float y_n, float t_n, float h, ODEFunction f, void* params)
{
    float k1 = h * f(t_n, y_n, params);
    float k2 = h * f(t_n + h * 0.5f, y_n + k1 * 0.5f, params);
    float k3 = h * f(t_n + h * 0.5f, y_n + k2 * 0.5f, params);
    float k4 = h * f(t_n + h, y_n + k3, params);

    return y_n + (1.0f / 6.0f) * (k1 + 2.0f * k2 + 2.0f * k3 + k4);
}


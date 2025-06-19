#include "simulation/solver/ODE45.h"

float ode45_solver(
    float y0, float t0, float t_end,
    float h_init, float tol,
    ODEFunction f, void* params)
{
    float t = t0;
    float y = y0;
    float h = h_init;

    const float safety = 0.9f;
    const float min_scale = 0.1f;
    const float max_scale = 5.0f;

    while (t < t_end) {
        if (t + h > t_end) {
            h = t_end - t;
        }

        // Dormand-Prince 5(4) coefficients (embedded)
        float k1 = h * f(t, y, params);
        float k2 = h * f(t + h * (1.0f / 5.0f), y + k1 * (1.0f / 5.0f), params);
        float k3 = h * f(t + h * (3.0f / 10.0f), y + k1 * (3.0f / 40.0f) + k2 * (9.0f / 40.0f), params);
        float k4 = h * f(t + h * (4.0f / 5.0f), y + k1 * (44.0f / 45.0f) - k2 * (56.0f / 15.0f) + k3 * (32.0f / 9.0f), params);
        float k5 = h * f(t + h * (8.0f / 9.0f), y + k1 * (19372.0f / 6561.0f) - k2 * (25360.0f / 2187.0f)
                                                  + k3 * (64448.0f / 6561.0f) - k4 * (212.0f / 729.0f), params);
        float k6 = h * f(t + h, y + k1 * (9017.0f / 3168.0f) - k2 * (355.0f / 33.0f)
                                      + k3 * (46732.0f / 5247.0f) + k4 * (49.0f / 176.0f)
                                      - k5 * (5103.0f / 18656.0f), params);
        float k7 = h * f(t + h, y + k1 * (35.0f / 384.0f) + 0.0f
                                      + k3 * (500.0f / 1113.0f) + k4 * (125.0f / 192.0f)
                                      - k5 * (2187.0f / 6784.0f) + k6 * (11.0f / 84.0f), params);

        // 5th order estimate
        float y5 = y + k1 * (35.0f / 384.0f) + k3 * (500.0f / 1113.0f)
                     + k4 * (125.0f / 192.0f) - k5 * (2187.0f / 6784.0f)
                     + k6 * (11.0f / 84.0f);

        // 4th order estimate
        float y4 = y + k1 * (5179.0f / 57600.0f) + k3 * (7571.0f / 16695.0f)
                     + k4 * (393.0f / 640.0f) - k5 * (92097.0f / 339200.0f)
                     + k6 * (187.0f / 2100.0f) + k7 * (1.0f / 40.0f);

        // Estimate local truncation error
        float err = fabsf(y5 - y4);

        // Accept step if error is within tolerance
        if (err <= tol || h < 1e-6f) {
            t += h;
            y = y5;
        } else {
        }

        // Update step size
        float scale = safety * powf(tol / (err + FLT_EPSILON), 0.2f);
        scale = fmaxf(min_scale, fminf(scale, max_scale));
        h *= scale;
    }

    return y;
}




#include "simulation/solver/solver.h"
#include "simulation/solver/runge_cutta.h"
#include "simulation/solver/ODE45.h"
#include "simulation/solver/RKF78.h"
#include <stdio.h>
#include <stdlib.h>


float solve(SolverType type, float tol, float y0, float t0, float t_end, float h_init, ODEFunction f, void* params)
{
    switch (type) {
    case SOLVER_RK4: {
        float t = t0;
        float y = y0;
        float h = h_init;
        while (t < t_end) {
            if (t + h > t_end)
                h = t_end - t;
            y = runge_kutta_4(y, t, h, f, params);
            t += h;
        }
        return y;
    }
    case SOLVER_ODE45:
        return ode45_solver(y0, t0, t_end, h_init, tol, f, params);
    case SOLVER_RKF78:
        return rkf78_solver(y0, t0, t_end, h_init, tol, f, params);
    default:
        fprintf(stderr, "[solve] Unknown solver type\n");
        exit(EXIT_FAILURE);
    }
}



float solve_step(SolverType type, float tol, float y_n, float t_n, float h, ODEFunction f, void* params) {
    switch (type) {
        case SOLVER_RK4:
            return runge_kutta_4(y_n, t_n, h, f, params);
        case SOLVER_ODE45:
            return ode45_solver(y_n, t_n, t_n + h, h, tol, f, params);
        case SOLVER_RKF78:
            return rkf78_solver(y_n, t_n, t_n + h, h, tol, f, params);
        default:
            fprintf(stderr, "[solve_step] Unknown solver type\n");
            exit(EXIT_FAILURE);
    }
}


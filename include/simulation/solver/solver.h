#ifndef SOLVER_H
#define SOLVER_H

#include <math.h>
#include <float.h>

typedef enum {
    SOLVER_RK4,
    SOLVER_ODE45,
    SOLVER_RKF78
} SolverType;



typedef float (*ODEFunction)(float t, float y, void* params);

float solve(SolverType type, float tol, float y0, float t0, float t_end, float h_init, ODEFunction f, void* params);
float solve_step(SolverType type, float tol, float y_n, float t_n, float h, ODEFunction f, void* params);


#endif // !SOLVER_H



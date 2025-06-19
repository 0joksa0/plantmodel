#ifndef ODE45_H
#define ODE45_H

#include "solver.h"

float ode45_solver(
    float y0, float t0, float t_end,
    float h_init, float tol,
    ODEFunction f, void* params);


#endif // !ODE45_H

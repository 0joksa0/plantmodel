#ifndef RKF78_H
#define RKF78_H

#include "solver.h"

float rkf78_solver(
    float y0, float t0, float t_end,
    float h_init, float tol,
    ODEFunction f, void* params);


#endif // !RKF78_H



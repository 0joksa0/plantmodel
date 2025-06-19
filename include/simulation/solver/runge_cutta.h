#ifndef RUNGE_CUTTA_H
#define RUNGE_CUTTA_H
#include "solver.h"

float runge_kutta_4(float y_n, float t_n, float h, ODEFunction f, void* params);


#endif

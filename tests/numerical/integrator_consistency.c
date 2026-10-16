#define SOLVER_REAL_AS double
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "solver.h"
static real_t f(real_t t, real_t y, void* u){ (void)t;(void)u; return -y; }
int gui_main_thread(void* a){(void)a;return 0;}
int main(void){
    SolverType kinds[3] = { SOLVER_RK4, SOLVER_ODE45, SOLVER_RKF78 };
    const char* nm[3] = {"RK4","ODE45","RKF78"};
    double tols[2] = {1e-4, 1e-6};
    for (int k=0;k<3;++k) for (int ti=0;ti<2;++ti){
        SolverRunConfig c; memset(&c,0,sizeof(c));
        c.type=kinds[k]; c.t0=0.0; c.t_end=1.0; c.h_init=0.01; c.tol=(real_t)tols[ti];
        real_t y = solver_scalar_run(&c, (real_t)1.0, f, NULL);
        double yy=(double)y;
        printf("%-6s tol=%.0e  y(1)=%.9f  exact=%.9f  relerr=%.3e  implied_rate=%.9f\n",
               nm[k], tols[ti], yy, exp(-1.0), fabs(yy-exp(-1.0))/exp(-1.0), -log(yy));
    }
    printf("1 - 82/840 = %.9f     1 - 41/840 = %.9f\n", 1.0-82.0/840.0, 1.0-41.0/840.0);
    return 0;
}

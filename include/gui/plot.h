#ifndef PLOT_H 
#define PLOT_H

#include <solver.h>
extern int total_steps;
extern real_t *sucrose;      /* niz[total_steps] */
extern real_t *starch;
extern real_t *ph;
extern real_t *partition;


void main_thread(void);
#endif 

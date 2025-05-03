#ifndef PLOT_H 
#define PLOT_H

extern int total_steps;
extern float *sucrose;      /* niz[total_steps] */
extern float *starch;
extern float *ph;
extern float *partition;


void main_thread(void);
#endif 

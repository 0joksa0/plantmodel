#ifndef PLOT_H
#define PLOT_H

#include <solver.h>

typedef struct {
    const real_t* sucrose;
    const real_t* starch;
    const real_t* photosynthesis;
    const real_t* total_biomass;
    const int* filled_steps;
    int total_steps;
    real_t dt_hours;
    real_t photoperiod;
    int days;
} GuiOutput;

typedef struct {
    const GuiOutput* output;
    const char* title;
    int target_fps;
    const char* summary_csv_path;
} GuiThreadArgs;

void* gui_main_thread(void* arg);

#endif

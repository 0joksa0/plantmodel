/* Quick λ_g sweep — 8h FJ, calibrated λ_sdi/λ_sb fixed, vary λ_g */
#include "model/input.h"
#include "model/model.h"
#include "simulation/simulation.h"
#include <stdio.h>
#include <string.h>

#define DAYS        30
#define DT_HOURS    REAL(0.1)

static const double LG[] = { 0.10, 0.20, 0.35, 0.50, 0.65 };
#define N_LG 5

int main(void)
{
    SimulationConfig cfg = simulation_default_config();
    cfg.days       = DAYS;
    cfg.dt_hours   = DT_HOURS;
    cfg.gui_enabled = 0;
    cfg.write_csv   = 0;

    for (int i = 0; i < N_LG; i++) {
        Input in = generate_input();
        in.photo.photo_mode                  = PHOTO_MODE_FARQUHAR_JARVIS;
        in.photo.photoperiod                 = REAL(8.0);
        in.photo.feedback_on_photosynthesis  = REAL(0.67);
        in.carbohydrates.lambda_sni          = REAL(0.13);
        in.growth.lambda_sb                  = REAL(0.005304);
        in.carbohydrates.lambda_sdi          = REAL(0.24980);
        in.carbohydrates.lambda_g            = (real_t)LG[i];

        SimulationResult res;
        memset(&res, 0, sizeof(res));
        simulate_days(&cfg, &in, &res);

        char fname[64];
        snprintf(fname, sizeof(fname), "/tmp/lg_%.2f.csv", LG[i]);
        FILE* f = fopen(fname, "w");
        fprintf(f, "zt_h,sucrose,starch\n");

        int spd   = (int)(REAL(24.0) / cfg.dt_hours);
        int start = res.filled_steps - spd;
        if (start < 0) start = 0;
        for (int j = start; j < res.filled_steps; j++) {
            double zt = (j - start) * (double)cfg.dt_hours;
            fprintf(f, "%.2f,%.6f,%.6f\n", zt, (double)res.sucrose[j], (double)res.starch[j]);
        }
        fclose(f);
        simulation_result_free(&res);
        printf("λ_g=%.2f  done → %s\n", LG[i], fname);
        fflush(stdout);
    }
    return 0;
}

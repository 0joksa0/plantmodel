/* Quick validation: simulate 8h FJ with vacuolar compartment, write last-day CSV */
#include "model/input.h"
#include "model/model.h"
#include "simulation/simulation.h"
#include <stdio.h>
#include <string.h>

int main(void)
{
    SimulationConfig cfg = simulation_default_config();
    cfg.days        = 30;
    cfg.dt_hours    = 0.1;
    cfg.gui_enabled = 0;
    cfg.write_csv   = 0;

    Input in = generate_input();
    in.photo.photo_mode                 = PHOTO_MODE_FARQUHAR_JARVIS;
    in.photo.photoperiod                = 8.0;
    in.photo.feedback_on_photosynthesis = 0.670;
    in.carbohydrates.lambda_sni         = 0.130;
    in.growth.lambda_sb                 = 0.005304;
    in.carbohydrates.lambda_sdi         = 0.24980;

    SimulationResult res;
    memset(&res, 0, sizeof(res));
    simulate_days(&cfg, &in, &res);

    int spd   = (int)(24.0 / cfg.dt_hours);
    int start = res.filled_steps - spd;
    if (start < 0) start = 0;

    puts("zt_h,sucrose_cyt,sucrose_vac,sucrose_total,starch");
    for (int j = start; j < res.filled_steps; j++) {
        double zt = (j - start) * (double)cfg.dt_hours;
        printf("%.2f,%.6f,%.6f,%.6f,%.6f\n",
               zt,
               (double)res.sucrose[j],
               (double)res.sucrose_vacuole[j],
               (double)res.sucrose_total[j],
               (double)res.starch[j]);
    }

    simulation_result_free(&res);
    return 0;
}

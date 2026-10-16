/* Generate last-day timeseries for all 4 photoperiods (FJ, calibrated params).
 * Writes results/ablation/timeseries/fj_Xh.csv with sucrose_total column. */
#include "model/input.h"
#include "model/model.h"
#include "simulation/simulation.h"
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>

typedef struct {
    int    pp;
    double lc;
    double lsni;
    double lsb;
    double lsdi;
    double k_vac_in;  /* calibrated k_vacuole_in per PP (k_out fixed at 0.08) */
} FJConfig;

/* k_vac_in calibrated by fj_sucrose_calibrate to hit reference mean sucrose_total
 * (Annunziata et al. 2017); k_vacuole_out = 0.08 h-1 fixed for all PP */
static const FJConfig CONFIGS[] = {
    {  4, 0.82, 0.16, 0.00384188, 0.05000000, 0.22994666 },
    {  6, 0.79, 0.15, 0.00448566, 0.30000000, 0.11925656 },
    {  8, 0.67, 0.13, 0.00530428, 0.24979858, 0.09649614 },
    { 12, 0.62, 0.08, 0.00591245, 0.05869751, 0.01134287 },
};

/* lambda_sb recalibrated against the same published RGR values under the
 * published growth conditions (docs/revision/light_dose_study.md).  The
 * production values above were bisected under a half-sine profile peaking at
 * 1200 umol m^-2 s^-1, whose daily light integral is 4.8x the one the plants
 * actually received, so they cannot be reused at the corrected irradiance. */
static const double LSB_CHAMBER[4] = { 0.0063494, 0.0069927, 0.0079587, 0.0082765 };

int main(int argc, char** argv)
{
    /* --chamber reproduces the growth conditions of Sulpice et al. 2014:
     * a square light profile at a constant 160 umol m^-2 s^-1 and a
     * 21 C light / 19 C dark thermoperiod.  Without the flag the tool behaves
     * exactly as before. */
    int chamber = 0;
    const char* out_dir = "results/ablation/timeseries";
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--chamber") == 0) {
            chamber = 1;
            out_dir = "results/ablation/timeseries_chamber";
        } else if (strcmp(argv[i], "--out-dir") == 0 && i + 1 < argc) {
            out_dir = argv[++i];
        } else {
            fprintf(stderr, "usage: %s [--chamber] [--out-dir DIR]\n", argv[0]);
            return 1;
        }
    }

    mkdir("results", 0755);
    mkdir("results/ablation", 0755);
    mkdir(out_dir, 0755);

    SimulationConfig cfg = simulation_default_config();
    cfg.days        = 30;
    cfg.dt_hours    = 0.1;
    cfg.gui_enabled = 0;
    cfg.write_csv   = 0;

    for (int i = 0; i < 4; i++) {
        const FJConfig* c = &CONFIGS[i];

        Input in = generate_input();
        in.photo.photo_mode                 = PHOTO_MODE_FARQUHAR_JARVIS;
        in.photo.photoperiod                = (double)c->pp;
        in.photo.feedback_on_photosynthesis = c->lc;
        in.carbohydrates.lambda_sni         = c->lsni;
        in.growth.lambda_sb                 = c->lsb;
        in.carbohydrates.lambda_sdi         = c->lsdi;
        in.carbohydrates.k_vacuole_in       = c->k_vac_in;
        in.carbohydrates.sucrose_vacuole    = (c->k_vac_in / 0.08) * in.carbohydrates.sucrose;

        if (chamber) {
            in.photo.light_profile                 = 1.0;   /* square wave */
            in.photo.max_PAR                       = 160.0; /* chamber irradiance */
            in.gas_exchange.leaf_temperature_day   = 21.0;
            in.gas_exchange.leaf_temperature_night = 19.0;
            in.gas_exchange.leaf_temperature       = 21.0;
            in.growth.lambda_sb                    = LSB_CHAMBER[i];
        }

        SimulationResult res;
        memset(&res, 0, sizeof(res));
        simulate_days(&cfg, &in, &res);

        char path[256];
        snprintf(path, sizeof(path), "%s/fj_%dh.csv", out_dir, c->pp);
        FILE* f = fopen(path, "w");
        if (!f) { perror(path); continue; }

        fprintf(f, "zt_h,sucrose_total,sucrose_cyt,sucrose_vac,starch\n");

        int spd   = (int)(24.0 / cfg.dt_hours);
        int start = res.filled_steps - spd;
        if (start < 0) start = 0;
        for (int j = start; j < res.filled_steps; j++) {
            double zt = (j - start) * (double)cfg.dt_hours;
            fprintf(f, "%.2f,%.6f,%.6f,%.6f,%.6f\n",
                    zt,
                    (double)res.sucrose_total[j],
                    (double)res.sucrose[j],
                    (double)res.sucrose_vacuole[j],
                    (double)res.starch[j]);
        }
        fclose(f);
        simulation_result_free(&res);
        printf("PP=%dh done → %s\n", c->pp, path);
        fflush(stdout);
    }
    return 0;
}

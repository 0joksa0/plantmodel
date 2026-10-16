/* 2D environment response surface: sweep CO2 × leaf_temperature for 4 PPs.
 * For each (CO2, T) point runs a 30-day simulation and outputs final RGR.
 * Output CSV: photoperiod,co2,temp,rgr
 */
#include "model/input.h"
#include "model/model.h"
#include "simulation/simulation.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

#define N_CO2   15
#define N_TEMP  15
#define N_DAYS  30

static const double CO2_MIN  = 200.0;
static const double CO2_MAX  = 900.0;
static const double TEMP_MIN = 8.0;
static const double TEMP_MAX = 38.0;

static const int    PPS[]   = { 4, 6, 8, 12 };
static const double LCS[]   = { 0.82, 0.79, 0.67, 0.62 };
static const double LSNIS[] = { 0.16, 0.15, 0.13, 0.08 };
static const double LSBS[]  = { 0.003832, 0.004480, 0.005301, 0.005867 };
static const double LSDIS[] = { 0.05,     0.30,     0.2507,   0.0605   };
static const double KINS[]  = { 0.22994666, 0.11925656, 0.09649614, 0.01134287 };

int main(void)
{
    SimulationConfig cfg = simulation_default_config();
    cfg.days        = N_DAYS;
    cfg.dt_hours    = 0.1;
    cfg.gui_enabled = 0;
    cfg.write_csv   = 0;

    printf("photoperiod,co2,temp,rgr\n");

    for (int pi = 0; pi < 4; pi++) {
        int    pp   = PPS[pi];
        double lc   = LCS[pi];
        double lsni = LSNIS[pi];
        double lsb  = LSBS[pi];
        double lsdi = LSDIS[pi];
        double kin  = KINS[pi];

        for (int ci = 0; ci < N_CO2; ci++) {
            double co2 = CO2_MIN + (CO2_MAX - CO2_MIN) * ci / (double)(N_CO2 - 1);

            for (int ti = 0; ti < N_TEMP; ti++) {
                double temp = TEMP_MIN + (TEMP_MAX - TEMP_MIN) * ti / (double)(N_TEMP - 1);

                Input in = generate_input();
                in.photo.photo_mode                 = PHOTO_MODE_FARQUHAR_JARVIS;
                in.photo.photoperiod                = (double)pp;
                in.photo.feedback_on_photosynthesis = lc;
                in.carbohydrates.lambda_sni         = lsni;
                in.growth.lambda_sb                 = lsb;
                in.carbohydrates.lambda_sdi         = lsdi;
                in.carbohydrates.k_vacuole_in       = kin;
                in.carbohydrates.sucrose_vacuole    = (kin / 0.08) * in.carbohydrates.sucrose;
                in.gas_exchange.ambient_CO2_concentration = co2;
                in.gas_exchange.leaf_temperature          = temp;

                double fw_start = in.growth.total_biomass;

                SimulationResult res;
                memset(&res, 0, sizeof(res));
                simulate_days(&cfg, &in, &res);
                simulation_result_free(&res);

                double fw_end = in.growth.total_biomass;
                double rgr = 0.0;
                if (fw_end > 0.0 && fw_start > 0.0)
                    rgr = (log(fw_end) - log(fw_start)) / (double)(N_DAYS - 1);

                printf("%d,%.2f,%.2f,%.6f\n", pp, co2, temp, rgr);
            }
        }
        fflush(stdout);
    }
    return 0;
}

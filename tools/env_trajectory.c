/* Dynamic environment simulation: CO2 and temperature evolve over N_DAYS.
 * Runs sequential 1-day simulations, passing final state as next initial state.
 * Outputs CSV: photoperiod,day,co2,temp,total_biomass,leaf_biomass,root_biomass,sucrose,starch
 */
#include "model/input.h"
#include "model/model.h"
#include "simulation/simulation.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

#define N_DAYS     90
#define N_WARMUP   20   /* warm-up days at baseline before environment shifts */

/* CO2 ramps from 280 (pre-industrial) to 750 ppm over the trajectory */
static double co2_schedule(int day)
{
    if (day < N_WARMUP) return 280.0;
    double t = (double)(day - N_WARMUP) / (double)(N_DAYS - N_WARMUP - 1);
    return 280.0 + 470.0 * t;  /* 280 → 750 ppm */
}

/* Temperature: baseline 22°C + seasonal sinusoidal swing ±8°C */
static double temp_schedule(int day)
{
    double phase = 2.0 * M_PI * day / (double)N_DAYS;
    return 22.0 + 8.0 * sin(phase - M_PI * 0.5);  /* starts cool, peaks mid-sim */
}

static const int PPS[]  = { 4, 6, 8, 12 };
static const double LCS[]  = { 0.82, 0.79, 0.67, 0.62 };
static const double LSNIS[] = { 0.16, 0.15, 0.13, 0.08 };
static const double LSBS[]  = { 0.003832, 0.004480, 0.005301, 0.005867 };
static const double LSDIS[] = { 0.05,     0.30,     0.2507,   0.0605   };
static const double KINS[]  = { 0.22994666, 0.11925656, 0.09649614, 0.01134287 };

int main(void)
{
    SimulationConfig cfg = simulation_default_config();
    cfg.days        = 1;
    cfg.dt_hours    = 0.1;
    cfg.gui_enabled = 0;
    cfg.write_csv   = 0;

    printf("photoperiod,day,co2,temp,total_biomass,leaf_biomass,root_biomass,sucrose,starch\n");

    for (int pi = 0; pi < 4; pi++) {
        int    pp   = PPS[pi];
        double lc   = LCS[pi];
        double lsni = LSNIS[pi];
        double lsb  = LSBS[pi];
        double lsdi = LSDIS[pi];
        double kin  = KINS[pi];

        Input in = generate_input();
        in.photo.photo_mode                 = PHOTO_MODE_FARQUHAR_JARVIS;
        in.photo.photoperiod                = (double)pp;
        in.photo.feedback_on_photosynthesis = lc;
        in.carbohydrates.lambda_sni         = lsni;
        in.growth.lambda_sb                 = lsb;
        in.carbohydrates.lambda_sdi         = lsdi;
        in.carbohydrates.k_vacuole_in       = kin;
        in.carbohydrates.sucrose_vacuole    = (kin / 0.08) * in.carbohydrates.sucrose;

        for (int day = 0; day < N_DAYS; day++) {
            double co2  = co2_schedule(day);
            double temp = temp_schedule(day);

            in.gas_exchange.ambient_CO2_concentration = co2;
            in.gas_exchange.leaf_temperature          = temp;

            SimulationResult res;
            memset(&res, 0, sizeof(res));
            simulate_days(&cfg, &in, &res);
            simulation_result_free(&res);

            printf("%d,%d,%.1f,%.2f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                   pp, day + 1,
                   co2, temp,
                   (double)in.growth.total_biomass,
                   (double)in.growth.leaf_biomass,
                   (double)in.growth.root_biomass,
                   (double)in.carbohydrates.sucrose,
                   (double)in.carbohydrates.starch);
        }
        fflush(stdout);
    }
    return 0;
}

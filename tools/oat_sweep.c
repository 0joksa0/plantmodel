/* One-at-a-time parameter sweep for model verification.
 * Each parameter swept over [lower, upper] (20 points) while others held nominal.
 * Output CSV: param,value,photoperiod,rgr,status
 *
 * `status` is one of:
 *   ok           -- trajectory finite, growth rate reported in `rgr`
 *   diverged     -- integration produced a non-finite state; `rgr` is empty
 *   nonpositive  -- finite trajectory but non-positive biomass; `rgr` is empty
 */
#include "model/input.h"
#include "model/model.h"
#include "simulation/simulation.h"
#include "util/fpclass.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

#define N_POINTS 20
#define N_DAYS   30
#define DT_HOURS 0.25   /* coarser step avoids solver stiffness at extreme params */

static const int    PPS[]   = { 4, 6, 8, 12 };
static const double LCS[]   = { 0.82, 0.79, 0.67, 0.62 };
static const double LSNIS[] = { 0.16, 0.15, 0.13, 0.08 };
static const double LSBS[]  = { 0.003832, 0.004480, 0.005301, 0.005867 };
static const double LSDIS[] = { 0.05,     0.30,     0.2507,   0.0605   };
static const double KINS[]  = { 0.22994666, 0.11925656, 0.09649614, 0.01134287 };

typedef struct {
    const char *name;
    double lower;
    double upper;
    int    log_scale;   /* 1 = log-uniform sweep */
} Param;

static const Param PARAMS[] = {
    { "max_nitrogen_uptake",          1.0,   15.0,  1 },
    { "max_phosphorus_uptake",        0.05,   1.5,  1 },
    { "leaf_water_potential",        -1.5,   -0.1,  0 },
    { "max_rubisco_carboxylation_rate", 20.0, 120.0, 1 },
    { "electron_transport_rate",      50.0,  200.0,  1 },
    { "respiration_frequency",         0.4,    1.2,  0 },
    { "leaf_temperature",             10.0,   35.0,  0 },
    { "ambient_CO2_concentration",   280.0,  800.0,  0 },
    { "relative_humidity",            30.0,   90.0,  0 },
    { "max_PAR",                     300.0, 2000.0,  0 },
};
#define N_PARAMS (int)(sizeof(PARAMS)/sizeof(PARAMS[0]))

static void set_param(Input *in, const char *name, double val)
{
    if (!strcmp(name, "max_nitrogen_uptake"))
        in->nutrients.max_nitrogen_uptake = val;
    else if (!strcmp(name, "max_phosphorus_uptake"))
        in->nutrients.max_phosphorus_uptake = val;
    else if (!strcmp(name, "leaf_water_potential"))
        in->gas_exchange.leaf_water_potential = val;
    else if (!strcmp(name, "max_rubisco_carboxylation_rate"))
        in->gas_exchange.max_rubisco_carboxylation_rate = val;
    else if (!strcmp(name, "electron_transport_rate"))
        in->gas_exchange.electron_transport_rate = val;
    else if (!strcmp(name, "respiration_frequency"))
        in->carbohydrates.respiration_frequency = val;
    else if (!strcmp(name, "leaf_temperature"))
        in->gas_exchange.leaf_temperature = val;
    else if (!strcmp(name, "ambient_CO2_concentration"))
        in->gas_exchange.ambient_CO2_concentration = val;
    else if (!strcmp(name, "relative_humidity"))
        in->gas_exchange.relative_humidity = val;
    else if (!strcmp(name, "max_PAR"))
        in->photo.max_PAR = val;
}

int main(void)
{
    SimulationConfig cfg = simulation_default_config();
    cfg.days        = N_DAYS;
    cfg.dt_hours    = 0.1;
    cfg.gui_enabled = 0;
    cfg.write_csv   = 0;
    cfg.solver_type = SOLVER_RK4;   /* fixed step: prevents adaptive step-size collapse at extreme params */

    printf("param,value,photoperiod,rgr,status\n");
    fflush(stdout);

    for (int pi = 0; pi < N_PARAMS; pi++) {
        const Param *p = &PARAMS[pi];
        fprintf(stderr, "param[%d] = %s\n", pi, p->name);

        for (int vi = 0; vi < N_POINTS; vi++) {
            double t = (double)vi / (double)(N_POINTS - 1);
            double val;
            if (p->log_scale)
                val = exp(log(p->lower) + t * (log(p->upper) - log(p->lower)));
            else
                val = p->lower + t * (p->upper - p->lower);

            for (int qi = 0; qi < 4; qi++) {
                int pp = PPS[qi];
                fprintf(stderr, "  vi=%d val=%.3f pp=%d\n", vi, val, pp); fflush(stderr);

                Input in = generate_input();
                in.photo.photo_mode                 = PHOTO_MODE_FARQUHAR_JARVIS;
                in.photo.photoperiod                = (double)pp;
                in.photo.feedback_on_photosynthesis = LCS[qi];
                in.carbohydrates.lambda_sni         = LSNIS[qi];
                in.growth.lambda_sb                 = LSBS[qi];
                in.carbohydrates.lambda_sdi         = LSDIS[qi];
                in.carbohydrates.k_vacuole_in       = KINS[qi];
                in.carbohydrates.sucrose_vacuole    = (KINS[qi] / 0.08) * in.carbohydrates.sucrose;

                set_param(&in, p->name, val);

                double fw_start = in.growth.total_biomass;
                SimulationResult res;
                memset(&res, 0, sizeof(res));
                simulate_days(&cfg, &in, &res);
                /* Read the divergence flag before the result buffers are freed. */
                int diverged = res.diverged;
                simulation_result_free(&res);

                double fw_end = in.growth.total_biomass;

                /* Classify the outcome explicitly.  A comparison against a NaN
                 * is false, so the historical `if (fw_end > 0.0)` guard silently
                 * mapped a diverged trajectory onto rgr = 0 -- indistinguishable
                 * from a plant that genuinely stopped growing.  Divergence is
                 * reported through `status` with an empty `rgr` field instead.
                 * Note that pm_isfinite() must be used rather than isfinite():
                 * the production build enables -ffast-math, under which the
                 * library macro is constant-folded to 1. */
                double      rgr    = 0.0;
                const char *status = "ok";
                if (diverged || !pm_isfinite(fw_end) || !pm_isfinite(fw_start)) {
                    status = "diverged";
                } else if (fw_end > 0.0 && fw_start > 0.0) {
                    rgr = (log(fw_end) - log(fw_start)) / (double)(N_DAYS - 1);
                } else {
                    status = "nonpositive";
                }

                if (strcmp(status, "ok") == 0)
                    printf("%s,%.6f,%d,%.6f,%s\n", p->name, val, pp, rgr, status);
                else
                    printf("%s,%.6f,%d,,%s\n", p->name, val, pp, status);
            }
        }
        fflush(stdout);
    }
    return 0;
}

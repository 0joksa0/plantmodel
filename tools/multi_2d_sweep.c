/* Sweep 4 biologically interesting parameter pairs × all 4 photoperiods.
 * Uses RK4 fixed step (prevents adaptive-solver collapse at extreme params).
 * Output CSV: pair,xname,yname,x,y,photoperiod,rgr,status
 *   status = ok          : finite trajectory, rgr is the reported growth rate
 *   status = diverged    : integration produced NaN/Inf; rgr field is empty
 *   status = nonpositive : finite but biomass did not stay positive; rgr = 0
 */
#include "model/input.h"
#include "model/model.h"
#include "simulation/simulation.h"
#include "util/fpclass.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

#define N_GRID  14
#define N_DAYS  30

static const int    PPS[]   = { 4, 6, 8, 12 };
static const double LCS[]   = { 0.82, 0.79, 0.67, 0.62 };
static const double LSNIS[] = { 0.16, 0.15, 0.13, 0.08 };
static const double LSBS[]  = { 0.003832, 0.004480, 0.005301, 0.005867 };
static const double LSDIS[] = { 0.05,     0.30,     0.2507,   0.0605   };
static const double KINS[]  = { 0.22994666, 0.11925656, 0.09649614, 0.01134287 };

typedef struct {
    const char *pair_name;
    const char *xname;  const char *yname;
    double xlo;  double xhi;  int xlog;
    double ylo;  double yhi;  int ylog;
} Pair;

static const Pair PAIRS[] = {
    /* 1. Nutrient co-limitation: N uptake × CO2 */
    { "N_CO2",
      "max_nitrogen_uptake",       "ambient_CO2_concentration",
       1.0,  15.0, 1,    280.0, 800.0, 0 },
    /* 2. Stomatal climate: temperature × relative humidity */
    { "T_RH",
      "leaf_temperature",          "relative_humidity",
      10.0,  35.0, 0,     30.0,  90.0, 0 },
    /* 3. Compound stress: water potential × temperature */
    { "psi_T",
      "leaf_water_potential",      "leaf_temperature",
      -1.5,  -0.1, 0,     10.0,  35.0, 0 },
    /* 4. Photosynthesis machinery: Vcmax × Jmax */
    { "Vcmax_Jmax",
      "max_rubisco_carboxylation_rate", "electron_transport_rate",
      20.0, 120.0, 1,     50.0, 200.0, 1 },
};
#define N_PAIRS (int)(sizeof(PAIRS)/sizeof(PAIRS[0]))

static void apply(Input *in, const char *name, double val)
{
    if      (!strcmp(name, "max_nitrogen_uptake"))
        in->nutrients.max_nitrogen_uptake = val;
    else if (!strcmp(name, "max_phosphorus_uptake"))
        in->nutrients.max_phosphorus_uptake = val;
    else if (!strcmp(name, "leaf_water_potential"))
        in->gas_exchange.leaf_water_potential = val;
    else if (!strcmp(name, "max_rubisco_carboxylation_rate"))
        in->gas_exchange.max_rubisco_carboxylation_rate = val;
    else if (!strcmp(name, "electron_transport_rate"))
        in->gas_exchange.electron_transport_rate = val;
    else if (!strcmp(name, "leaf_temperature"))
        in->gas_exchange.leaf_temperature = val;
    else if (!strcmp(name, "ambient_CO2_concentration"))
        in->gas_exchange.ambient_CO2_concentration = val;
    else if (!strcmp(name, "relative_humidity"))
        in->gas_exchange.relative_humidity = val;
    else if (!strcmp(name, "max_PAR"))
        in->photo.max_PAR = val;
}

static double sample(double lo, double hi, int lg, int i)
{
    double t = (double)i / (double)(N_GRID - 1);
    return lg ? exp(log(lo) + t * (log(hi) - log(lo)))
              : lo + t * (hi - lo);
}

int main(void)
{
    SimulationConfig cfg = simulation_default_config();
    cfg.days        = N_DAYS;
    cfg.dt_hours    = 0.1;
    cfg.gui_enabled = 0;
    cfg.write_csv   = 0;
    cfg.solver_type = SOLVER_RK4;

    printf("pair,xname,yname,x,y,photoperiod,rgr,status\n");
    fflush(stdout);

    for (int pi = 0; pi < N_PAIRS; pi++) {
        const Pair *p = &PAIRS[pi];
        fprintf(stderr, "pair[%d] = %s\n", pi, p->pair_name);

        for (int xi = 0; xi < N_GRID; xi++) {
            double xval = sample(p->xlo, p->xhi, p->xlog, xi);
            for (int yi = 0; yi < N_GRID; yi++) {
                double yval = sample(p->ylo, p->yhi, p->ylog, yi);

                for (int qi = 0; qi < 4; qi++) {
                    int pp = PPS[qi];

                    Input in = generate_input();
                    in.photo.photo_mode                 = PHOTO_MODE_FARQUHAR_JARVIS;
                    in.photo.photoperiod                = (double)pp;
                    in.photo.feedback_on_photosynthesis = LCS[qi];
                    in.carbohydrates.lambda_sni         = LSNIS[qi];
                    in.growth.lambda_sb                 = LSBS[qi];
                    in.carbohydrates.lambda_sdi         = LSDIS[qi];
                    in.carbohydrates.k_vacuole_in       = KINS[qi];
                    in.carbohydrates.sucrose_vacuole    = (KINS[qi]/0.08)*in.carbohydrates.sucrose;

                    apply(&in, p->xname, xval);
                    apply(&in, p->yname, yval);

                    double fw0 = in.growth.total_biomass;
                    SimulationResult res;
                    memset(&res, 0, sizeof(res));
                    simulate_days(&cfg, &in, &res);
                    int diverged = res.diverged;
                    simulation_result_free(&res);

                    /*
                     * Outcome classification. The previous version collapsed
                     * every non-growing outcome to rgr = 0.0 via a bare
                     * `if (total_biomass > 0.0)` guard. That is wrong for
                     * divergent integrations: at extreme (Vcmax, Jmax) or
                     * (N_s, C_a) combinations the plant outgrows the soil
                     * nutrient supply, nitrogen crosses zero, and the
                     * affinity ODEs (n_nmax = max_N/(N + max_N)) hit their
                     * pole, so the trajectory becomes NaN. `NaN > 0.0` is
                     * false, so those cells were written as rgr = 0 and were
                     * indistinguishable from genuine zero growth — they
                     * appeared on response surfaces as a plateau of "high
                     * capacity kills the plant", which the model does not say.
                     *
                     * Non-finite outcomes are now reported as such (empty rgr
                     * field, status = diverged) so downstream plots can mask
                     * them instead of drawing them as data.
                     */
                    double Bt = in.growth.total_biomass;
                    /* Two independent checks: the in-run guard record, and a
                     * fold-proof test of the reported biomass. pm_isfinite,
                     * not isfinite — -ffast-math folds the library macro to a
                     * constant. See include/util/fpclass.h */
                    int fin = !diverged
                              && pm_isfinite(Bt) && pm_isfinite(fw0);
                    const char *status;
                    double rgr = 0.0;

                    if (!fin) {
                        status = "diverged";
                    } else if (Bt > 0.0 && fw0 > 0.0) {
                        rgr = (log(Bt) - log(fw0)) / (N_DAYS - 1);
                        status = "ok";
                    } else {
                        status = "nonpositive";
                    }

                    if (fin)
                        printf("%s,%s,%s,%.6f,%.6f,%d,%.6f,%s\n",
                               p->pair_name, p->xname, p->yname,
                               xval, yval, pp, rgr, status);
                    else
                        printf("%s,%s,%s,%.6f,%.6f,%d,,%s\n",
                               p->pair_name, p->xname, p->yname,
                               xval, yval, pp, status);
                }
            }
        }
        fflush(stdout);
    }
    return 0;
}

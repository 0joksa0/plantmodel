/* nsweep -- nitrogen saturation curve against the measured series, over the
 * range the experiment actually used, for a grid of (I_N_max, k_n).
 *
 * The measured series (Sulpice et al. 2014, N_data / nitrogen_response_reference)
 * gives tissue N at hydroponic supplies of 3, 25, 100, 400 and 2000 uM. The
 * model's soil term is umol N cm^-3 = umol/mL, so 1 uM = 1e-3 umol cm^-3 and the
 * measured range is 0.003-2.0. The shipped k_n = 0.125 puts half-saturation at
 * 125 uM; the measurement's midpoint sits near 400 uM.
 *
 * Scored on SHAPE, normalised to the lowest supply, so a difference in absolute
 * pool size does not enter -- the model's N pool and a tissue percentage are not
 * the same quantity and must not be compared unnormalised.
 */
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "model/model.h"
#include "config/config.h"
#include "simulation/simulation.h"

#define NOBS 5
static const double SUP_UM[NOBS]  = { 3, 25, 100, 400, 2000 };
static const double TISSUE[NOBS]  = { 1.920, 2.160, 3.231667, 4.868333, 4.985 };

static int run_one(double soil, double imax, double kn, int chamber, double pp, double *pool)
{
    SimulationConfig cfg; Input in;
    cfg = simulation_default_config();
    in  = generate_input();
    if (config_load_file("config/default.config", &cfg, &in) != 0) return 1;
    cfg.days = 30; cfg.write_csv = 0;
    in.photo.photo_mode = (int)PHOTO_MODE_FARQUHAR_JARVIS;
    in.photo.photoperiod = (real_t)pp;
    in.photo.feedback_on_photosynthesis = (real_t)(pp<=4.5?0.82: pp<=6.5?0.79: pp<=8.5?0.67:0.62);
    in.carbohydrates.lambda_sni = (real_t)(pp<=4.5?0.16: pp<=6.5?0.15: pp<=8.5?0.13:0.08);
    /* lambda_sb = 1/eta, the reciprocal of the construction cost.  It is no
     * longer a fitted value: see growth_cost() in src/model/carbohydrates.c. */
    in.growth.lambda_sb = REAL(0.00154419);
    in.nutrients.nitrogen_soil_content = (real_t)soil;
    in.nutrients.max_nitrogen_uptake = (real_t)imax;
    in.nutrients.Michaelis_Menten_constant_nitrogen = (real_t)kn;
    if (chamber) {
        in.photo.max_PAR = REAL(160.0);
        in.gas_exchange.leaf_temperature = REAL(21.0);
#ifdef PLANTMODEL_HAS_DAYNIGHT_TEMP
        in.gas_exchange.leaf_temperature_day = REAL(21.0);
        in.gas_exchange.leaf_temperature_night = REAL(19.0);
#endif
#ifdef PLANTMODEL_HAS_LIGHT_PROFILE
        in.photo.light_profile = REAL(1.0);
#endif
    }
    SimulationResult r; memset(&r, 0, sizeof(r));
    simulate_days(&cfg, &in, &r);
    int bad = (r.filled_steps != r.total_steps);
    *pool = bad ? NAN : (double)in.nutrients.nitrogen;
    if (!isfinite(*pool)) bad = 1;
    simulation_result_free(&r);
    return bad;
}

int main(int argc, char** argv)
{
    int chamber = 1;
    const char* out = "out/nsweep.csv";
    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "--production")) chamber = 0;
        else if (!strcmp(argv[i], "--out") && i + 1 < argc) out = argv[++i];
    }
    FILE* fp = fopen(out, "w");
    fprintf(fp, "imax,kn,chamber,photoperiod,supply_uM,pool,pool_norm,tissue_norm,ok\\n");

    /* Adopted pair only, swept across all four photoperiods. The grid over
     * (Imax, kn) was measured earlier and is not repeated here; kn = 0.125 is
     * T3 (Tedone Eq. 2.2.22) and is not a free parameter. */
    double IMAX[] = { 6.44 };
    double KN[]   = { 0.125 };
    double PPS[]  = { 4.0, 6.0, 8.0, 12.0 };

    for (unsigned a = 0; a < sizeof IMAX / sizeof *IMAX; ++a)
    for (unsigned b = 0; b < sizeof KN / sizeof *KN; ++b)
    for (unsigned q = 0; q < sizeof PPS / sizeof *PPS; ++q) {
        double pool[NOBS]; int allok = 1;
        for (int k = 0; k < NOBS; ++k) {
            if (run_one(SUP_UM[k] / 1000.0, IMAX[a], KN[b], chamber, PPS[q], &pool[k]))
                { allok = 0; pool[k] = NAN; }
        }
        for (int k = 0; k < NOBS; ++k)
            fprintf(fp, "%.4f,%.4f,%d,%g,%g,%.6f,%.6f,%.6f,%d\n",
                    IMAX[a], KN[b], chamber, PPS[q], SUP_UM[k], pool[k],
                    allok ? pool[k] / pool[0] : NAN,
                    TISSUE[k] / TISSUE[0], allok);
        fflush(fp);
        printf("pp=%4.1f Imax=%5.2f kn=%.3f  %s\n", PPS[q], IMAX[a], KN[b], allok ? "ok" : "FAILED");
    }
    fclose(fp);
    return 0;
}

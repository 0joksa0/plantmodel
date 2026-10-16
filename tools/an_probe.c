/* an_probe -- the area-based assimilation rate the conversion factor multiplies.
 *
 * WHY.  The Sulpice route to the conversion factor is a ratio: the measured
 * per-gram rate divided by the model's own area-based rate under the SAME light.
 * The denominator has been quoted as An = 5.154 umol CO2 m-2 s-1 "under chamber
 * conditions", meaning a SQUARE 160 umol m-2 s-1.  But at this baseline commit
 * update_light_conditions() imposes a HALF-SINE (photosynthesis.c:341), and the
 * square-wave profile field does not exist here, so --chamber gives a half-sine
 * peaking at 160 whose mean is 160*2/pi = 101.9.  The denominator therefore
 * depends on which light schedule is meant, and so does the factor.
 *
 * This tool evaluates the model's own photosynthesis routine at a FIXED PAR, so
 * the denominator is measured rather than assumed.  It reports the raw area-based
 * An by dividing out the conversion factor that farquhar_photosynthesis() has
 * already applied: An_area = photosynthesis / kappa.
 */

#include "config/config.h"
#include "model/input.h"
#include "model/model.h"
#include "simulation/simulation.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

static void ensure_dir(const char* p)
{
    char b[512]; snprintf(b, sizeof(b), "%s", p);
    for (char* q = b + 1; *q; ++q) if (*q == '/') { *q = 0; mkdir(b, 0755); *q = '/'; }
}

int main(int argc, char** argv)
{
    const char* out = "results/attribution/an_probe.csv";
    double temp = 21.0;
    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "--out") && i + 1 < argc)  out = argv[++i];
        else if (!strcmp(argv[i], "--temp") && i + 1 < argc) temp = atof(argv[++i]);
    }
    ensure_dir(out);
    FILE* fp = fopen(out, "w");
    if (!fp) return 1;
    fprintf(fp, "PAR,leaf_temp,kappa,An_area_umol_CO2_m2_s,A_mass_umol_hex_gFW_h\n");

    const double PARS[] = { 100.0, 101.86, 120.0, 160.0, 200.0, 250.0, 300.0,
                            500.0, 764.0, 1200.0 };
    const int NP = (int)(sizeof(PARS) / sizeof(PARS[0]));

    for (int j = 0; j < NP; ++j) {
        SimulationConfig cfg = simulation_default_config();
        Input in = generate_input();
        if (config_load_file("config/default.config", &cfg, &in) != 0) return 1;
        in.photo.photo_mode = (int)PHOTO_MODE_FARQUHAR_JARVIS;
        in.gas_exchange.leaf_temperature = (real_t)temp;
        in.photo.light = REAL(1.0);
        in.photo.light_PAR = (real_t)PARS[j];

        /* farquhar_photosynthesis() overwrites light_PAR only via
         * update_light_conditions(), which we deliberately do not call: the
         * point is to hold PAR fixed. */
        real_t A = farquhar_photosynthesis(&in);
        double kappa = 1.10 * (double)in.gas_exchange.specific_leaf_area / 0.025;
        fprintf(fp, "%.2f,%.1f,%.4f,%.6f,%.6f\n",
                PARS[j], temp, kappa, (double)A / kappa, (double)A);
    }
    fclose(fp);
    printf("wrote %s\n", out);
    return 0;
}

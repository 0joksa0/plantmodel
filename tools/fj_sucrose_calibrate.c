/* Calibrate k_vacuole_in per photoperiod to match reference mean sucrose_total.
 * Bisects k_in in [lo, hi] until mean(sucrose_total over last diel cycle) = target.
 * All other params (λ_sb, λ_sdi, λ_sni, lc) held at existing FJ ablation values.
 *
 * Reference sucrose means from Annunziata et al. 2017 (µmol/gFW):
 *   4h: 2.4075   6h: 2.3420   8h: 2.7826   12h: 2.5861
 */
#include "model/input.h"
#include "model/model.h"
#include "simulation/simulation.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

/* Reference means (ZT 0-22, every 2h, 12 points each) */
static const double REF_SUC_MEAN[4] = { 2.4075, 2.3420, 2.7826, 2.5861 };

typedef struct {
    int    pp;
    double lc;
    double lsni;
    double lsb;
    double lsdi;
} FJCfg;

static const FJCfg CFGS[4] = {
    {  4, 0.82, 0.16, 0.00384188, 0.05000000 },
    {  6, 0.79, 0.15, 0.00448566, 0.30000000 },
    {  8, 0.67, 0.13, 0.00530428, 0.24979858 },
    { 12, 0.62, 0.08, 0.00591245, 0.05869751 },
};

static double sim_mean_sucrose(const FJCfg* c, double k_in, double k_out)
{
    SimulationConfig cfg = simulation_default_config();
    cfg.days        = 40;
    cfg.dt_hours    = 0.1;
    cfg.gui_enabled = 0;
    cfg.write_csv   = 0;

    Input in = generate_input();
    in.photo.photo_mode                 = PHOTO_MODE_FARQUHAR_JARVIS;
    in.photo.photoperiod                = (double)c->pp;
    in.photo.feedback_on_photosynthesis = c->lc;
    in.carbohydrates.lambda_sni         = c->lsni;
    in.growth.lambda_sb                 = c->lsb;
    in.carbohydrates.lambda_sdi         = c->lsdi;
    in.carbohydrates.k_vacuole_in       = k_in;
    in.carbohydrates.k_vacuole_out      = k_out;
    /* initialize vacuolar pool at algebraic SS for fast convergence */
    in.carbohydrates.sucrose_vacuole    = (k_in / k_out) * in.carbohydrates.sucrose;

    SimulationResult res;
    memset(&res, 0, sizeof(res));
    simulate_days(&cfg, &in, &res);

    int spd   = (int)(24.0 / cfg.dt_hours);
    int start = res.filled_steps - spd;
    if (start < 0) start = 0;

    double sum = 0.0;
    int    n   = 0;
    for (int j = start; j < res.filled_steps; j++) {
        sum += (double)res.sucrose_total[j];
        n++;
    }
    simulation_result_free(&res);
    return (n > 0) ? sum / n : 0.0;
}

int main(void)
{
    const double K_OUT = 0.08;
    const double TOL   = 5e-4;  /* µmol/gFW */

    printf("%-6s  %-12s  %-12s  %-12s\n", "PP", "k_in", "mean_sucrose", "target");
    printf("------  ------------  ------------  ------------\n");
    fflush(stdout);

    for (int i = 0; i < 4; i++) {
        const FJCfg* c      = &CFGS[i];
        double        target = REF_SUC_MEAN[i];

        double lo = 1e-3, hi = 4.0;
        double k_in = 0.0;

        double vlo = sim_mean_sucrose(c, lo, K_OUT);
        double vhi = sim_mean_sucrose(c, hi, K_OUT);

        if (vlo >= target) {
            k_in = lo;
        } else if (vhi <= target) {
            k_in = hi;
        } else {
            for (int it = 0; it < 50; it++) {
                k_in     = 0.5 * (lo + hi);
                double v = sim_mean_sucrose(c, k_in, K_OUT);
                if (fabs(v - target) < TOL) break;
                if (v < target) lo = k_in;
                else            hi = k_in;
            }
        }

        double final_v = sim_mean_sucrose(c, k_in, K_OUT);
        printf("%-6dh  %-12.8f  %-12.4f  %-12.4f\n",
               c->pp, k_in, final_v, target);
        fflush(stdout);
    }

    return 0;
}

/* Compute RGR for all 4 photoperiods using calibrated FJ + vacuolar k_in params.
 * Prints CSV: pp,lit_rgr,model_rgr,rel_err_pct */
#include "model/input.h"
#include "model/model.h"
#include "simulation/simulation.h"
#include <stdio.h>
#include <string.h>

static const int DAYS = 30;

typedef struct {
    int    pp;
    double lit_rgr;
    double lc;
    double lsni;
    double lsb;
    double lsdi;
    double k_vac_in;
} FJCfg;

static const FJCfg CFGS[] = {
    {  4, 0.0680, 0.82, 0.16, 0.00384188, 0.05000000, 0.22994666 },
    {  6, 0.1135, 0.79, 0.15, 0.00448566, 0.30000000, 0.11925656 },
    {  8, 0.1708, 0.67, 0.13, 0.00530428, 0.24979858, 0.09649614 },
    { 12, 0.2600, 0.62, 0.08, 0.00591245, 0.05869751, 0.01134287 },
};

int main(void)
{
    SimulationConfig cfg = simulation_default_config();
    cfg.days        = DAYS;
    cfg.dt_hours    = 0.1;
    cfg.gui_enabled = 0;
    cfg.write_csv   = 0;

    printf("pp,lit_rgr,model_rgr,rel_err_pct\n");

    for (int i = 0; i < 4; i++) {
        const FJCfg* c = &CFGS[i];

        Input in = generate_input();
        in.photo.photo_mode                 = PHOTO_MODE_FARQUHAR_JARVIS;
        in.photo.photoperiod                = (double)c->pp;
        in.photo.feedback_on_photosynthesis = c->lc;
        in.carbohydrates.lambda_sni         = c->lsni;
        in.growth.lambda_sb                 = c->lsb;
        in.carbohydrates.lambda_sdi         = c->lsdi;
        in.carbohydrates.k_vacuole_in       = c->k_vac_in;
        in.carbohydrates.sucrose_vacuole    = (c->k_vac_in / 0.08) * in.carbohydrates.sucrose;

        double fw_start = in.growth.total_biomass;

        SimulationResult res;
        memset(&res, 0, sizeof(res));
        simulate_days(&cfg, &in, &res);
        simulation_result_free(&res);

        double fw_end   = in.growth.total_biomass;
        double rgr      = compute_rgr(fw_start, fw_end, (double)(DAYS - 1));
        double rel_err  = (rgr - c->lit_rgr) / c->lit_rgr * 100.0;

        printf("%d,%.6f,%.6f,%.2f\n", c->pp, c->lit_rgr, rgr, rel_err);
        fflush(stdout);
    }
    return 0;
}

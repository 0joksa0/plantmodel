/* Sweep a SINGLE lambda_sb across all 4 photoperiods (chamber conditions),
 * compute RGR from biomass, report error vs literature targets.
 * Purpose: test whether one global lambda_sb reproduces the RGR-vs-PP
 * relationship without per-photoperiod refitting. */
#include "model/input.h"
#include "model/model.h"
#include "simulation/simulation.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

typedef struct { int pp; double lc; double lsni; double lsdi; double k_vac_in; double lit_rgr; } Cfg;
/* Table A.1 fixed (lc, lsni) + calibrated lsdi/k_vac_in from fj_diel_runner;
 * lambda_sb is the swept single value. lit_rgr = published target. */
static const Cfg C[4] = {
  {  4, 0.82, 0.16, 0.05000000, 0.22994666, 0.0680 },
  {  6, 0.79, 0.15, 0.30000000, 0.11925656, 0.1135 },
  {  8, 0.67, 0.13, 0.24979858, 0.09649614, 0.1708 },
  { 12, 0.62, 0.08, 0.05869751, 0.01134287, 0.2600 },
};

static double run_rgr(const Cfg* c, double lsb, int days) {
    SimulationConfig cfg = simulation_default_config();
    cfg.days=days; cfg.dt_hours=0.1; cfg.gui_enabled=0; cfg.write_csv=0;
    Input in = generate_input();
    in.photo.photo_mode                 = PHOTO_MODE_FARQUHAR_JARVIS;
    in.photo.photoperiod                = (double)c->pp;
    in.photo.feedback_on_photosynthesis = c->lc;
    in.carbohydrates.lambda_sni         = c->lsni;
    in.growth.lambda_sb                 = lsb;
    in.carbohydrates.lambda_sdi         = c->lsdi;
    in.carbohydrates.k_vacuole_in       = c->k_vac_in;
    in.carbohydrates.sucrose_vacuole    = (c->k_vac_in/0.08)*in.carbohydrates.sucrose;
    /* chamber conditions (the experiment's own irradiance) */
    in.photo.light_profile                 = 1.0;
    in.photo.max_PAR                       = 160.0;
    in.gas_exchange.leaf_temperature_day   = 21.0;
    in.gas_exchange.leaf_temperature_night = 19.0;
    in.gas_exchange.leaf_temperature       = 21.0;
    double fw0 = in.growth.total_biomass;
    SimulationResult res; memset(&res,0,sizeof(res));
    simulate_days(&cfg,&in,&res);
    simulation_result_free(&res);
    double fw1 = in.growth.total_biomass;
    if (!(fw1>0.0) || !isfinite(fw1)) return NAN;
    return compute_rgr(fw0, fw1, (double)(days-1));
}

int main(int argc, char** argv) {
    int days=30;
    double lo=0.0040, hi=0.0110; int N=71;
    printf("lambda_sb");
    for (int i=0;i<4;i++) printf(",rgr_%dh,err_%dh", C[i].pp, C[i].pp);
    printf(",MAPE_pct,nfinite\n");
    double best_lsb=NAN, best_mape=1e30;
    for (int g=0; g<N; ++g) {
        double lsb = lo + (hi-lo)*g/(double)(N-1);
        double sae=0; int nf=0, nok=0;
        double rgr[4], err[4];
        for (int i=0;i<4;i++){
            rgr[i]=run_rgr(&C[i], lsb, days);
            if (isfinite(rgr[i])){ err[i]=100.0*fabs(rgr[i]-C[i].lit_rgr)/C[i].lit_rgr; sae+=err[i]; nok++; nf++; }
            else err[i]=NAN;
        }
        double mape = nok? sae/nok : NAN;
        printf("%.6f", lsb);
        for (int i=0;i<4;i++) printf(",%.6f,%.3f", rgr[i], err[i]);
        printf(",%.4f,%d\n", mape, nf);
        if (nf==4 && mape<best_mape){ best_mape=mape; best_lsb=lsb; }
    }
    fprintf(stderr, "BEST single lambda_sb=%.6f  MAPE=%.3f%% (all 4 finite)\n", best_lsb, best_mape);
    return 0;
}

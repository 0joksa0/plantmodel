#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "model/model.h"
#include "config/config.h"
#include "simulation/simulation.h"
int gui_main_thread(void* a){(void)a;return 0;}
int main(int argc,char**argv){
    SimulationConfig cfg = simulation_default_config();
    Input in = generate_input();
    if (config_load_file("config/default.config",&cfg,&in)!=0){printf("cfg fail\n");return 1;}
    double pp=(argc>1)?atof(argv[1]):12.0;
    in.photo.photoperiod=pp; cfg.days=35;
    SimulationResult r; simulate_days(&cfg,&in,&r);
    if(r.diverged){printf("pp=%.0f DIVERGED\n",pp);return 0;}
    double eff_max_starch = (double)max_starch(in.carbohydrates.max_starch_degradation_rate, in.photo.photoperiod);
    double cap_suc = (double)in.carbohydrates.max_sucrose;
    long nsuc=0, nst=0; double smax=0, stmax=0, fmin=1e9;
    for(int i=0;i<r.filled_steps;i++){
        double su=r.sucrose[i], st=r.starch[i];
        if(su>smax) smax=su; if(st>stmax) stmax=st;
        if(su>cap_suc) nsuc++;
        if(st>eff_max_starch) nst++;
        double f = 0.62 + (1.0-0.62)*((eff_max_starch-st)/eff_max_starch);
        if(f<fmin) fmin=f;
    }
    printf("pp=%2.0f | sucrose: max=%.4f cap=%.2f  ABOVE_CAP=%5ld/%d (%.1f%%)"
           " | starch: max=%.3f eff_cap=%.1f above=%ld | min f_starch=%.4f\n",
           pp, smax, cap_suc, nsuc, r.filled_steps, 100.0*nsuc/r.filled_steps,
           stmax, eff_max_starch, nst, fmin);
    simulation_result_free(&r);
    return 0;
}

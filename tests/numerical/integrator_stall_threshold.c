#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>
#include <sys/wait.h>
#include "config/config.h"
#include "model/model.h"
#include "simulation/simulation.h"
int gui_main_thread(void* a){(void)a;return 0;}
static SimulationConfig g_cfg; static Input g_base;
static void onalarm(int s){(void)s; _exit(3);}
/* child-process probe: returns 0 = completed, 3 = stalled */
static int try_run(double pp,double lsb,double lsdi,int days,int timeout_s,double*rgr){
    fflush(stdout);
    pid_t k=fork();
    if(k==0){
        signal(SIGALRM,onalarm); alarm(timeout_s);
        g_cfg.days=days;
        Input in=g_base;
        in.photo.photo_mode=(int)PHOTO_MODE_FARQUHAR_JARVIS;
        in.photo.photoperiod=(real_t)pp;
        in.carbohydrates.lambda_sdi=(real_t)lsdi;
        in.growth.lambda_sb=(real_t)lsb;
        real_t fw0=in.growth.total_biomass;
        SimulationResult r; memset(&r,0,sizeof(r));
        simulate_days(&g_cfg,&in,&r);
        double g=(double)compute_rgr(fw0,in.growth.total_biomass,(real_t)(days-1));
        FILE*f=fopen("/dev/shm/stall_rgr","w"); fprintf(f,"%.10g\n",g); fclose(f);
        _exit(0);
    }
    int st; waitpid(k,&st,0);
    int code=WIFEXITED(st)?WEXITSTATUS(st):-1;
    *rgr=NAN;
    if(code==0){FILE*f=fopen("/dev/shm/stall_rgr","r"); if(f){fscanf(f,"%lf",rgr);fclose(f);} }
    return code;
}
int main(void){
    g_cfg=simulation_default_config(); g_base=generate_input();
    config_load_file("config/default.config",&g_cfg,&g_base);
    g_cfg.gui_enabled=0; g_cfg.write_csv=0;
    double pps[4]={4,6,8,12};
    printf("Bisecting the lambda_sdi stall threshold (lambda_sb = 0.0066631, days = 8, 15 s cap)\n\n");
    for(int i=0;i<4;i++){
        double lo=0.02, hi=1.00, rgr;
        /* confirm lo completes and hi stalls */
        int clo=try_run(pps[i],0.0066630853,lo,8,15,&rgr);
        int chi=try_run(pps[i],0.0066630853,hi,8,15,&rgr);
        printf("pp=%2.0f  lo(%.2f)=%s  hi(%.2f)=%s",pps[i],lo,clo==0?"ok":"STALL",hi,chi==0?"ok":"STALL");
        if(clo!=0){printf("   -> stalls across the whole range\n");continue;}
        if(chi==0){printf("   -> no stall in range\n");continue;}
        for(int it=0;it<12;it++){
            double mid=0.5*(lo+hi);
            if(try_run(pps[i],0.0066630853,mid,8,15,&rgr)==0) lo=mid; else hi=mid;
        }
        printf("   threshold in [%.4f, %.4f]\n",lo,hi);
        fflush(stdout);
    }
    return 0;
}

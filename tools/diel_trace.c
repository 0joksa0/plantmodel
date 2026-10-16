/* Diel trajectory export for the standard report panel.
 *
 * Writes one CSV per photoperiod with the columns the reporting layer expects:
 * t_h, starch, sucrose, sucrose_vac, photosynthesis -- sampled over the final
 * 24 h of a 30-day run, so the plant is on its limit cycle. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <model/model.h>
#include <config/config.h>
#include <simulation/simulation.h>
#include <simulation/dynamics.h>
static double PP[8]={4,6,8,12};
static int NPP_RUN=4;
int main(int argc,char**argv){
  const char* dir="out"; const char* tag="trace";
  double lsb=0.00154419; int chamber=0, days=30;
  /* Negative means "keep the per-photoperiod literature value" (the historical
   * behaviour).  A non-negative value applies the SAME number at every
   * photoperiod, which is what a configuration with no per-photoperiod fitting
   * requires -- and what model_report's --shared-lc/--shared-lsni already do.
   * Without these flags this tool silently overwrote the config value, so the
   * carbohydrate panels described a different model from the RGR table. */
  double shared_lc=-1.0, shared_lsni=-1.0;
  for(int i=1;i<argc;i++){
    if(!strcmp(argv[i],"--dir")&&i+1<argc) dir=argv[++i];
    else if(!strcmp(argv[i],"--tag")&&i+1<argc) tag=argv[++i];
    else if(!strcmp(argv[i],"--lambda-sb")&&i+1<argc) lsb=atof(argv[++i]);
    else if(!strcmp(argv[i],"--chamber")) chamber=1;
    else if(!strcmp(argv[i],"--days")&&i+1<argc) days=atoi(argv[++i]);
    else if(!strcmp(argv[i],"--shared-lc")&&i+1<argc) shared_lc=atof(argv[++i]);
    else if(!strcmp(argv[i],"--shared-lsni")&&i+1<argc) shared_lsni=atof(argv[++i]);
    else if(!strcmp(argv[i],"--pp")&&i+1<argc){ PP[0]=atof(argv[++i]); NPP_RUN=1; }
  }
  mkdir(dir,0755);
  for(int i=0;i<NPP_RUN;i++){
    SimulationConfig c=simulation_default_config(); Input in=generate_input();
    if(config_load_file("config/default.config",&c,&in)!=0){
      fprintf(stderr,"cannot read config\n"); return 1; }
    c.gui_enabled=0; c.write_csv=0; c.days=days;
    in.photo.photo_mode=(int)PHOTO_MODE_FARQUHAR_JARVIS;
    in.photo.photoperiod=(real_t)PP[i];
    /* The per-photoperiod literature values that every other runner sets.
     * Omitting them here silently produced a different model: starch peaked at
     * 1800 instead of 82 umol C6 gFW-1. */
    in.photo.feedback_on_photosynthesis = shared_lc >= 0.0 ? (real_t)shared_lc :
        (real_t)(PP[i]<=4.5?0.82: PP[i]<=6.5?0.79: PP[i]<=8.5?0.67:0.62);
    in.carbohydrates.lambda_sni = shared_lsni >= 0.0 ? (real_t)shared_lsni :
        (real_t)(PP[i]<=4.5?0.16: PP[i]<=6.5?0.15: PP[i]<=8.5?0.13:0.08);
    in.growth.lambda_sb=(real_t)lsb;
    if(chamber){
      in.photo.light_profile=REAL(1.0); in.photo.max_PAR=REAL(160.0);
      in.gas_exchange.leaf_temperature_day=REAL(21.0);
      in.gas_exchange.leaf_temperature_night=REAL(19.0);
    }
    SimulationResult r={0}; simulate_days(&c,&in,&r);
    char path[512]; snprintf(path,sizeof path,"%s/%s_%dh.csv",dir,tag,(int)PP[i]);
    FILE* f=fopen(path,"w");
    if(!f){ fprintf(stderr,"cannot write %s\n",path); return 1; }
    fprintf(f,"t_h,starch,sucrose,sucrose_vac,photosynthesis\n");
    int spd=(int)lround(24.0/(double)c.dt_hours);
    int c0=r.filled_steps-spd; if(c0<0) c0=0;
    for(int k=c0;k<r.filled_steps;k++)
      fprintf(f,"%.4f,%.6f,%.6f,%.6f,%.6f\n",(k-c0)*(double)c.dt_hours,
        (double)r.starch[k],(double)r.sucrose[k],
        (double)r.sucrose_vacuole[k],(double)r.ph[k]);
    fclose(f);
    printf("%s  filled=%d/%d %s\n",path,r.filled_steps,r.total_steps,
      r.filled_steps==r.total_steps?"OK":"DIVERGED");
    simulation_result_free(&r);
  }
  return 0;
}

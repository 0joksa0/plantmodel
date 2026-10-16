/* Shared-parameter ceiling and held-out prediction for the recalibrated model.
 *
 * PROTOCOL (identical to the four-row comparison table in
 * docs/revision/respiration_and_recalibration.md, so the rows are comparable):
 *
 *   - ONE shared parameter set for all photoperiods: lambda_sb, lambda_sdr,
 *     lambda_sdi, lambda_sni.  Photoperiod is the only thing that differs
 *     between the four runs.  No per-photoperiod freedom anywhere.
 *   - Nelder-Mead, three deterministic starts, best of three reported.
 *   - Objective: sum of relative squared errors in RGR, dusk starch and dusk
 *     sucrose_total, over the photoperiods in the training set.
 *   - Chamber configuration: light_profile = 1 (square), max_PAR = 160,
 *     leaf_temperature_day = 21, leaf_temperature_night = 19.
 *   - kappa = 3.25 as specific_leaf_area = 0.025 * (3.25/1.10).
 *   - 30 days; dusk sampled at ZT = photoperiod on the final day.
 *
 * lambda_c (feedback_on_photosynthesis) is held at its per-photoperiod Table
 * A.1 values, as in every previous harness, and is NOT fitted here.  It is a
 * photosynthesis feedback and it is inert in FARQUHAR_JARVIS mode -- the audit
 * measured gross A_n identical to seven significant figures across the whole
 * scenario set -- so it costs no effective degrees of freedom.  Reported for
 * completeness rather than used.
 *
 * MODES
 *   --mode ceiling   fit on all four photoperiods; write the ceiling row.
 *   --mode heldout   every informative split: four 3-train/1-test and six
 *                    2-train/2-test.  Fit on train, FREEZE, predict test.
 *   --mode diel      dump the day-30 diel trace at the fitted parameters, for
 *                    the respiration-phasing figure.
 *
 * --release-clamp raises max_sucrose out of the operating range at run time.
 * It exists so that this same tool, compiled against the pre-change tree, can
 * reproduce the "+ both" predecessor row under exactly this protocol.  Against
 * the current tree the RMIN is gone from night_efficiency_starch() and the
 * flag only perturbs the Michaelis denominators, so it must be left off.
 *
 * COMPILE FLAGS: needs BOTH -DSOLVER_REAL_AS=double AND
 * -DSOLVER_REAL_IS_DOUBLE.  With only the first, RMAX/RMIN resolve to the
 * single-precision fmaxf/fminf and this is not the shipped model.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "model/model.h"
#include "config/config.h"
#include "simulation/simulation.h"
#include "simulation/dynamics.h"

/* No gui_main_thread stub in the CMake build: that target is built from
 * ${SRC}, which already includes src/gui/plot.c, so a stub would be a
 * duplicate definition.  Standalone builds that exclude src/gui -- notably the
 * before/after comparison, which compiles this same file against the pristine
 * pre-change tree under resp-build/before/ -- do need one, and ask for it with
 * -DPLANTMODEL_STANDALONE_GUI_STUB. */
#ifdef PLANTMODEL_STANDALONE_GUI_STUB
int gui_main_thread(void* a){(void)a;return 0;}
#endif

#define KAPPA 3.25
#define NPP   4
#define NPAR  5

static const int    PP[NPP]      = {4, 6, 8, 12};
static const double LC[NPP]      = {0.82, 0.79, 0.67, 0.62};
/* Sulpice et al. 2014: RGR from the growth analysis, dusk starch and dusk
 * sucrose_total from mmc4 "Means" at ZT = photoperiod. */
static const double RGR_OBS[NPP] = {0.0680, 0.1135, 0.1708, 0.2600};
static const double ST_OBS[NPP]  = {29.6749, 54.2372, 56.8778, 58.7337};
static const double SU_OBS[NPP]  = {1.722, 3.710, 3.418, 2.344};

static int g_release_clamp = 0;

typedef struct { double rgr, st, su; int ok; } Out;

static Out sim(int i, const double* p)
{
    Out o = {0,0,0,0};
    SimulationConfig cfg = simulation_default_config();
    Input in = generate_input();
    if (config_load_file("config/default.config", &cfg, &in) != 0) return o;
    cfg.days = 30;
    cfg.gui_enabled = 0;
    cfg.write_csv = 0;
    in.photo.photo_mode  = (int)PHOTO_MODE_FARQUHAR_JARVIS;
    in.photo.photoperiod = (real_t)PP[i];
    in.photo.feedback_on_photosynthesis = (real_t)LC[i];
    in.carbohydrates.lambda_sdr = (real_t)p[1];
    in.carbohydrates.lambda_sdi = (real_t)p[2];
    in.carbohydrates.lambda_sni = (real_t)p[3];
    in.growth.lambda_sb = (real_t)p[0];
    in.carbohydrates.sucrose_loading_frequency = (real_t)p[4];
    in.gas_exchange.specific_leaf_area = (real_t)(0.025 * (KAPPA / 1.10));
    in.photo.light_profile = REAL(1.0);
    in.photo.max_PAR = REAL(160.0);
    in.gas_exchange.leaf_temperature_day   = REAL(21.0);
    in.gas_exchange.leaf_temperature_night = REAL(19.0);
    in.gas_exchange.leaf_temperature       = REAL(21.0);
    if (g_release_clamp) in.carbohydrates.max_sucrose = REAL(1000.0);

    real_t fw0 = in.growth.total_biomass;
    SimulationResult r;
    memset(&r, 0, sizeof(r));
    simulate_days(&cfg, &in, &r);
    if (r.diverged) { simulation_result_free(&r); return o; }
    o.rgr = (double)compute_rgr(fw0, in.growth.total_biomass, (real_t)(cfg.days-1));
    int spd = (int)(24.0 / cfg.dt_hours);
    int idx = (cfg.days-1)*spd + (int)(PP[i]/cfg.dt_hours);
    if (idx >= r.filled_steps) idx = r.filled_steps - 1;
    o.st = (double)r.starch[idx];
    o.su = (double)r.sucrose_total[idx];
    o.ok = isfinite(o.rgr) && isfinite(o.st) && isfinite(o.su)
           && o.st < 1e5 && o.st > -1e-6 && o.rgr > 0.0;
    simulation_result_free(&r);
    return o;
}

/* Bit mask over photoperiods: bit i selects PP[i]. */
static int g_train_mask = 0x0F;

static int in_bounds(const double* p)
{
    return  p[0] >= 1e-5 && p[0] <= 0.05
        &&  p[1] >= 0.05 && p[1] <= 12.0
        &&  p[2] >= 0.005 && p[2] <= 0.60
        &&  p[3] >= 0.005 && p[3] <= 0.60
        /* Export coefficient.  The bounds are deliberately wide: the model
         * imposes its own, much tighter, admissibility window on this
         * parameter (measured with probe/window.c as slf in [6.3, 7.2] at the
         * ceiling lambdas -- bounded below by the 12 h photoperiod, whose
         * sucrose pool has no steady state if export is too weak, and above by
         * the 4 h photoperiod, which drains if export is too strong).  Letting
         * the optimiser see the whole range and be stopped by divergence keeps
         * that window a measured property of the model rather than something
         * written into the harness. */
        &&  p[4] >= 0.20  && p[4] <= 20.0;
}

/* Mean relative squared error over the selected photoperiods, three targets
 * each.  Mean rather than sum so that a 2-photoperiod split is on the same
 * scale as a 4-photoperiod one and train/test SSE can be compared directly. */
static double score_mask(const Out* o, int mask, int* nseen)
{
    double s = 0.0; int n = 0;
    for (int i = 0; i < NPP; ++i) {
        if (!((mask >> i) & 1)) continue;
        if (!o[i].ok) { if (nseen) *nseen = 0; return 1e6; }
        double er = (o[i].rgr - RGR_OBS[i]) / RGR_OBS[i];
        double es = (o[i].st  - ST_OBS[i])  / ST_OBS[i];
        double eu = (o[i].su  - SU_OBS[i])  / SU_OBS[i];
        s += er*er + es*es + eu*eu; ++n;
    }
    if (nseen) *nseen = n;
    return n ? s / n : 1e6;
}

static void run_all(const double* p, Out* o)
{
    #pragma omp parallel for schedule(dynamic) num_threads(4)
    for (int i = 0; i < NPP; ++i) o[i] = sim(i, p);
}

static double obj(const double* p)
{
    if (!in_bounds(p)) return 1e6;
    Out o[NPP];
    /* Only the training photoperiods are simulated; the test ones must not
     * influence the fit in any way, not even through a shared failure flag. */
    #pragma omp parallel for schedule(dynamic) num_threads(4)
    for (int i = 0; i < NPP; ++i) {
        if ((g_train_mask >> i) & 1) o[i] = sim(i, p);
        else { o[i].ok = 0; o[i].rgr = o[i].st = o[i].su = 0.0; }
    }
    return score_mask(o, g_train_mask, NULL);
}

static void nm(double* x, int iters)
{
    const int n = NPAR, m = NPAR + 1;
    double S[NPAR+1][NPAR], F[NPAR+1];
    const double step[NPAR] = {0.0010, 0.9, 0.06, 0.05, 0.8};
    for (int i = 0; i < m; ++i) {
        memcpy(S[i], x, n*sizeof(double));
        if (i > 0) S[i][i-1] += step[i-1];
        F[i] = obj(S[i]);
    }
    for (int it = 0; it < iters; ++it) {
        int hi=0, lo=0;
        for (int i=1;i<m;++i){ if(F[i]>F[hi])hi=i; if(F[i]<F[lo])lo=i; }
        double c[NPAR]={0,0,0,0};
        for (int i=0;i<m;++i) if(i!=hi) for(int j=0;j<n;++j) c[j]+=S[i][j]/(m-1);
        double t[NPAR], ft;
        for (int j=0;j<n;++j) t[j]=c[j]+(c[j]-S[hi][j]);
        ft = obj(t);
        if (ft < F[lo]) {
            double e[NPAR], fe;
            for (int j=0;j<n;++j) e[j]=c[j]+2.0*(c[j]-S[hi][j]);
            fe = obj(e);
            if (fe < ft) { memcpy(t,e,sizeof(t)); ft = fe; }
        } else if (ft >= F[hi]) {
            for (int j=0;j<n;++j) t[j]=c[j]+0.5*(S[hi][j]-c[j]);
            ft = obj(t);
        }
        if (ft < F[hi]) { memcpy(S[hi],t,sizeof(t)); F[hi]=ft; }
        else for (int i=0;i<m;++i) if(i!=lo) {
            for (int j=0;j<n;++j) S[i][j]=0.5*(S[i][j]+S[lo][j]);
            F[i]=obj(S[i]);
        }
    }
    int lo=0; for(int i=1;i<m;++i) if(F[i]<F[lo]) lo=i;
    memcpy(x, S[lo], n*sizeof(double));
}

/* Three deterministic starts, spread across the FEASIBLE region so that a bad
 * basin is distinguishable from a real ceiling.
 *
 * These are not arbitrary.  Once the respiratory term is applied to the total
 * sucrose pool the drain is 2.25x larger, and the export coefficient that
 * balances it moves with it: the whole neighbourhood of the old starts (slf
 * around 2) now fails to integrate at 12 h, so every simplex launched there
 * dies on the divergence sentinel before it can move.  probe/feas.c scans a
 * coarse grid for combinations admissible at all four photoperiods (112 of
 * 5760) and these three bracket the feasible set it found, with the first
 * being its best point.  Starting inside the feasible region is a property of
 * the harness, not of the model: the reported optimum is still whatever the
 * simplex converges to, and the three starts still have to agree. */
static const double STARTS[3][NPAR] = {
    {0.0030, 0.35, 0.15, 0.005, 6.00},
    {0.0026, 0.61, 0.31, 0.005, 6.50},
    {0.0036, 0.15, 0.05, 0.015, 5.50},
};

static double fit(int mask, double* best, int iters, int verbose)
{
    g_train_mask = mask;
    double bf = 1e9;
    for (int s = 0; s < 3; ++s) {
        double x[NPAR]; memcpy(x, STARTS[s], sizeof(x));
        nm(x, iters);
        double f = obj(x);
        if (verbose)
            fprintf(stderr, "   start %d -> lsb=%.7f lsdr=%.4f lsdi=%.4f lsni=%.4f slf=%.4f  J=%.5f\n",
                    s, x[0], x[1], x[2], x[3], x[4], f);
        if (f < bf) { bf = f; memcpy(best, x, NPAR*sizeof(double)); }
    }
    return bf;
}

static void mape(const Out* o, int mask, double* mr, double* ms, double* mu)
{
    double a=0,b=0,c=0; int n=0;
    for (int i=0;i<NPP;++i) {
        if (!((mask>>i)&1)) continue;
        a += fabs(100.0*(o[i].rgr-RGR_OBS[i])/RGR_OBS[i]);
        b += fabs(100.0*(o[i].st -ST_OBS[i]) /ST_OBS[i]);
        c += fabs(100.0*(o[i].su -SU_OBS[i]) /SU_OBS[i]);
        ++n;
    }
    *mr = n?a/n:NAN; *ms = n?b/n:NAN; *mu = n?c/n:NAN;
}

static void mask_label(int mask, char* buf)
{
    buf[0] = 0; int first = 1;
    for (int i=0;i<NPP;++i) if ((mask>>i)&1) {
        char t[8]; snprintf(t,sizeof t,"%s%d", first?"":"+", PP[i]);
        strcat(buf,t); first = 0;
    }
    if (first) strcpy(buf,"-");
}

static void mode_ceiling(const char* outdir, int iters)
{
    double best[NPAR];
    fprintf(stderr, "[ceiling] fitting all four photoperiods\n");
    double f = fit(0x0F, best, iters, 1);
    Out o[NPP]; g_train_mask = 0x0F; run_all(best, o);
    double mr, ms, mu; mape(o, 0x0F, &mr, &ms, &mu);
    /* SSE on the brief's convention: summed, not averaged, over 4 photoperiods */
    double sse = f * 4.0;

    printf("\n=== CEILING, one shared parameter set, kappa=%.2f, chamber ===\n", KAPPA);
    printf("lambda_sb=%.7f lambda_sdr=%.5f lambda_sdi=%.5f lambda_sni=%.5f slf=%.5f\n",
           best[0], best[1], best[2], best[3], best[4]);
    printf("RGR MAPE %.1f%%   starch MAPE %.1f%%   sucrose MAPE %.1f%%   SSE %.3f\n\n",
           mr, ms, mu, sse);
    printf("%4s %10s %10s %8s | %9s %9s %8s | %8s %8s %8s\n",
           "pp","RGR","obs","err%","starch","obs","err%","sucrose","obs","err%");
    for (int i=0;i<NPP;++i)
        printf("%3dh %10.6f %10.6f %8.2f | %9.2f %9.2f %8.1f | %8.2f %8.2f %8.1f\n",
               PP[i], o[i].rgr, RGR_OBS[i], 100.0*(o[i].rgr-RGR_OBS[i])/RGR_OBS[i],
               o[i].st, ST_OBS[i], 100.0*(o[i].st-ST_OBS[i])/ST_OBS[i],
               o[i].su, SU_OBS[i], 100.0*(o[i].su-SU_OBS[i])/SU_OBS[i]);

    char path[512];
    snprintf(path, sizeof path, "%s/ceiling_fit.csv", outdir);
    FILE* fp = fopen(path, "w");
    if (fp) {
        fprintf(fp, "photoperiod,rgr_model,rgr_obs,rgr_err_pct,starch_model,starch_obs,"
                    "starch_err_pct,sucrose_model,sucrose_obs,sucrose_err_pct\n");
        for (int i=0;i<NPP;++i)
            fprintf(fp, "%d,%.8f,%.8f,%.4f,%.6f,%.6f,%.4f,%.6f,%.6f,%.4f\n",
                    PP[i], o[i].rgr, RGR_OBS[i], 100.0*(o[i].rgr-RGR_OBS[i])/RGR_OBS[i],
                    o[i].st, ST_OBS[i], 100.0*(o[i].st-ST_OBS[i])/ST_OBS[i],
                    o[i].su, SU_OBS[i], 100.0*(o[i].su-SU_OBS[i])/SU_OBS[i]);
        fclose(fp);
    }
    snprintf(path, sizeof path, "%s/ceiling_params.csv", outdir);
    fp = fopen(path, "w");
    if (fp) {
        fprintf(fp, "configuration,lambda_sb,lambda_sdr,lambda_sdi,lambda_sni,sucrose_loading_frequency,"
                    "rgr_mape_pct,starch_mape_pct,sucrose_mape_pct,sse\n");
        fprintf(fp, "%s,%.8f,%.6f,%.6f,%.6f,%.6f,%.3f,%.3f,%.3f,%.4f\n",
                g_release_clamp ? "clamp_released_runtime" : "current_source",
                best[0], best[1], best[2], best[3], best[4], mr, ms, mu, sse);
        fclose(fp);
    }
}

static void mode_heldout(const char* outdir, int iters)
{
    /* Every informative split: |train| = 3 (4 splits) and |train| = 2 (6). */
    int masks[10]; int nm_ = 0;
    for (int m = 1; m < 16; ++m) {
        int bits = __builtin_popcount(m);
        if (bits == 3 || bits == 2) masks[nm_++] = m;
    }

    char path[512];
    snprintf(path, sizeof path, "%s/heldout_matrix.csv", outdir);
    FILE* fp = fopen(path, "w");
    if (fp) fprintf(fp, "train,test,n_train,n_test,lambda_sb,lambda_sdr,lambda_sdi,lambda_sni,sucrose_loading_frequency,"
                        "train_J,test_J,train_rgr_mape,train_starch_mape,train_sucrose_mape,"
                        "test_rgr_mape,test_starch_mape,test_sucrose_mape\n");
    snprintf(path, sizeof path, "%s/heldout_per_photoperiod.csv", outdir);
    FILE* fq = fopen(path, "w");
    if (fq) fprintf(fq, "train,test,photoperiod,role,rgr_model,rgr_obs,starch_model,starch_obs,"
                        "sucrose_model,sucrose_obs\n");

    printf("\n=== HELD-OUT PREDICTION ===\n");
    printf("%-10s %-10s %9s %9s | %8s %8s %8s | %8s %8s %8s\n",
           "train","test","train_J","test_J",
           "trRGR%","trStar%","trSuc%","teRGR%","teStar%","teSuc%");

    for (int k = 0; k < nm_; ++k) {
        int tr = masks[k], te = (~tr) & 0x0F;
        char ltr[32], lte[32]; mask_label(tr, ltr); mask_label(te, lte);
        fprintf(stderr, "[heldout] train %s test %s\n", ltr, lte);
        double best[NPAR];
        double trJ = fit(tr, best, iters, 0);
        /* FREEZE, then predict.  Nothing below touches the parameters. */
        Out o[NPP]; g_train_mask = 0x0F; run_all(best, o);
        int nte = 0;
        double teJ = score_mask(o, te, &nte);
        double tmr,tms,tmu,emr,ems,emu;
        mape(o, tr, &tmr,&tms,&tmu);
        mape(o, te, &emr,&ems,&emu);
        printf("%-10s %-10s %9.4f %9.4f | %8.1f %8.1f %8.1f | %8.1f %8.1f %8.1f\n",
               ltr, lte, trJ, teJ, tmr,tms,tmu, emr,ems,emu);
        fflush(stdout);
        if (fp) fprintf(fp, "%s,%s,%d,%d,%.8f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,"
                            "%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n",
                        ltr, lte, __builtin_popcount(tr), __builtin_popcount(te),
                        best[0],best[1],best[2],best[3],best[4], trJ, teJ,
                        tmr,tms,tmu, emr,ems,emu);
        if (fq) for (int i=0;i<NPP;++i)
            fprintf(fq, "%s,%s,%d,%s,%.8f,%.8f,%.6f,%.6f,%.6f,%.6f\n",
                    ltr, lte, PP[i], ((tr>>i)&1)?"train":"test",
                    o[i].rgr, RGR_OBS[i], o[i].st, ST_OBS[i], o[i].su, SU_OBS[i]);
    }
    if (fp) fclose(fp);
    if (fq) fclose(fq);
}

/* Day-30 diel trace of the respiratory flux, for the phasing figure.
 * maintenance respiration is recomputed here from the sampled cytosolic
 * sucrose and the light state, which is exactly what the RHS does. */
/* Evaluate one frozen parameter set at every photoperiod, reporting the
 * divergence flag explicitly.
 *
 * The held-out matrix reports a sentinel objective when a test photoperiod
 * fails to integrate, which is correct for the fit but hides WHY.  A failed
 * prediction can mean the frozen parameters are genuinely inadmissible at that
 * photoperiod, or merely that the run hit the step budget.  Those are
 * different findings and the difference has to be measurable, so this mode
 * prints the divergence kind alongside the values. */
static void mode_predict(const double* p)
{
    printf("\n=== FROZEN PARAMETER SET, evaluated at every photoperiod ===\n");
    printf("lambda_sb=%.7f lambda_sdr=%.6f lambda_sdi=%.6f lambda_sni=%.6f slf=%.6f\n",
           p[0], p[1], p[2], p[3], p[4]);
    printf("%4s %10s %10s %10s %12s %10s %10s\n",
           "pp","rgr","starch","sucrose","diverged","stalled","left_dom");
    for (int i = 0; i < NPP; ++i) {
        SimulationConfig cfg = simulation_default_config();
        Input in = generate_input();
        if (config_load_file("config/default.config", &cfg, &in) != 0) continue;
        cfg.days = 30; cfg.gui_enabled = 0; cfg.write_csv = 0;
        in.photo.photo_mode  = (int)PHOTO_MODE_FARQUHAR_JARVIS;
        in.photo.photoperiod = (real_t)PP[i];
        in.photo.feedback_on_photosynthesis = (real_t)LC[i];
        in.growth.lambda_sb = (real_t)p[0];
        in.carbohydrates.lambda_sdr = (real_t)p[1];
        in.carbohydrates.lambda_sdi = (real_t)p[2];
        in.carbohydrates.lambda_sni = (real_t)p[3];
        in.carbohydrates.sucrose_loading_frequency = (real_t)p[4];
        in.gas_exchange.specific_leaf_area = (real_t)(0.025*(KAPPA/1.10));
        in.photo.light_profile = REAL(1.0);
        in.photo.max_PAR = REAL(160.0);
        in.gas_exchange.leaf_temperature_day   = REAL(21.0);
        in.gas_exchange.leaf_temperature_night = REAL(19.0);
        in.gas_exchange.leaf_temperature       = REAL(21.0);
        if (g_release_clamp) in.carbohydrates.max_sucrose = REAL(1000.0);
        real_t fw0 = in.growth.total_biomass;
        SimulationResult r; memset(&r,0,sizeof(r));
        simulate_days(&cfg, &in, &r);
        int spd = (int)(24.0/cfg.dt_hours);
        int idx = (cfg.days-1)*spd + (int)(PP[i]/cfg.dt_hours);
        if (idx >= r.filled_steps) idx = r.filled_steps>0 ? r.filled_steps-1 : 0;
        double rgr = (double)compute_rgr(fw0, in.growth.total_biomass,
                                         (real_t)(cfg.days-1));
        printf("%3dh %10.6f %10.3f %10.3f %12d %10d %10d\n",
               PP[i], rgr, (double)r.starch[idx], (double)r.sucrose_total[idx],
               r.diverged, r.stalled, r.left_domain);
        simulation_result_free(&r);
    }
}

static void mode_diel(const char* outdir, const double* p)
{
    char path[512];
    snprintf(path, sizeof path, "%s/diel_respiration.csv", outdir);
    FILE* fp = fopen(path, "w");
    if (!fp) return;
    fprintf(fp, "photoperiod,ZT,light,sucrose_cyt,sucrose_total,starch,"
                "r_m_times_S,gated_resp,Rd_light_equiv\n");
    for (int i = 0; i < NPP; ++i) {
        SimulationConfig cfg = simulation_default_config();
        Input in = generate_input();
        if (config_load_file("config/default.config", &cfg, &in) != 0) continue;
        cfg.days = 30; cfg.gui_enabled = 0; cfg.write_csv = 0;
        in.photo.photo_mode  = (int)PHOTO_MODE_FARQUHAR_JARVIS;
        in.photo.photoperiod = (real_t)PP[i];
        in.photo.feedback_on_photosynthesis = (real_t)LC[i];
        in.carbohydrates.lambda_sdr = (real_t)p[1];
        in.carbohydrates.lambda_sdi = (real_t)p[2];
        in.carbohydrates.lambda_sni = (real_t)p[3];
        in.growth.lambda_sb = (real_t)p[0];
        in.carbohydrates.sucrose_loading_frequency = (real_t)p[4];
        in.gas_exchange.specific_leaf_area = (real_t)(0.025 * (KAPPA / 1.10));
        in.photo.light_profile = REAL(1.0);
        in.photo.max_PAR = REAL(160.0);
        in.gas_exchange.leaf_temperature_day   = REAL(21.0);
        in.gas_exchange.leaf_temperature_night = REAL(19.0);
        in.gas_exchange.leaf_temperature       = REAL(21.0);
        if (g_release_clamp) in.carbohydrates.max_sucrose = REAL(1000.0);
        double rm = (double)in.carbohydrates.respiration_frequency;
        /* Rd on the per-gram basis: the FvCB dark respiration, scaled by the
         * same conversion the assimilation rate goes through. */
        double Rd = (double)in.gas_exchange.respiration_rate * KAPPA;

        SimulationResult r; memset(&r, 0, sizeof(r));
        simulate_days(&cfg, &in, &r);
        if (r.diverged) { simulation_result_free(&r); continue; }
        int spd = (int)(24.0 / cfg.dt_hours);
        int base = (cfg.days-1)*spd;
        for (int k = 0; k < spd && base+k < r.filled_steps; ++k) {
            double zt = k * (double)cfg.dt_hours;
            int light = (zt < (double)PP[i]);
            /* The respiratory substrate is the TOTAL sucrose pool -- that is
             * the quantity r_m was measured against (whole-tissue extract),
             * and it is what maintenance_respiration() multiplies.  Reporting
             * r_m * S_cyt here would understate the flux by the vacuolar
             * factor of 2.25 and disagree with the model it is documenting. */
            double S    = (double)r.sucrose[base+k];
            double Stot = (double)r.sucrose_total[base+k];
            fprintf(fp, "%d,%.4f,%d,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f\n",
                    PP[i], zt, light, S, Stot,
                    (double)r.starch[base+k],
                    rm*Stot, (light?0.0:1.0)*rm*Stot, light?Rd:0.0);
        }
        simulation_result_free(&r);
    }
    fclose(fp);
    fprintf(stderr, "[diel] written\n");
}

/* 24 h integrals of every term in the sucrose budget on day 30, plus the
 * carbon actually incorporated into new biomass over the same day.  This is
 * the diagnostic that decides whether a residual is a respiration problem or
 * a growth-accounting problem: if the terms balance but the pools are still
 * too full, the sink that is missing is not respiratory. */
static void mode_budget(const char* outdir, const double* p)
{
    char path[512];
    snprintf(path, sizeof path, "%s/carbon_budget.csv", outdir);
    FILE* fp = fopen(path, "w");
    if (!fp) return;
    fprintf(fp, "photoperiod,A_to_sucrose,starch_deg,uptake_cost,transport_cost,"
                "maint_resp,export_load,vac_net,dS_pool,residual,"
                "biomass_C_incorporated,Cdens_implied,rgr\n");
    printf("\n=== 24 h SUCROSE BUDGET, day 30 [umol hexose gFW-1 d-1] ===\n");
    printf("%4s %10s %9s %8s %8s %9s %10s %8s | %9s %10s %10s\n",
           "pp","A->suc","starchdeg","uptake","transp","maint","export","vac",
           "dPool","resid","biomassC");
    for (int i = 0; i < NPP; ++i) {
        SimulationConfig cfg = simulation_default_config();
        Input in = generate_input();
        if (config_load_file("config/default.config", &cfg, &in) != 0) continue;
        cfg.days = 30; cfg.gui_enabled = 0; cfg.write_csv = 0;
        in.photo.photo_mode  = (int)PHOTO_MODE_FARQUHAR_JARVIS;
        in.photo.photoperiod = (real_t)PP[i];
        in.photo.feedback_on_photosynthesis = (real_t)LC[i];
        in.carbohydrates.lambda_sdr = (real_t)p[1];
        in.carbohydrates.lambda_sdi = (real_t)p[2];
        in.carbohydrates.lambda_sni = (real_t)p[3];
        in.growth.lambda_sb = (real_t)p[0];
        in.carbohydrates.sucrose_loading_frequency = (real_t)p[4];
        in.gas_exchange.specific_leaf_area = (real_t)(0.025*(KAPPA/1.10));
        in.photo.light_profile = REAL(1.0);
        in.photo.max_PAR = REAL(160.0);
        in.gas_exchange.leaf_temperature_day   = REAL(21.0);
        in.gas_exchange.leaf_temperature_night = REAL(19.0);
        in.gas_exchange.leaf_temperature       = REAL(21.0);
        if (g_release_clamp) in.carbohydrates.max_sucrose = REAL(1000.0);
        real_t fw0 = in.growth.total_biomass;
        SimulationResult r; memset(&r,0,sizeof(r));
        simulate_days(&cfg, &in, &r);
        if (r.diverged) { simulation_result_free(&r); continue; }
        double rgr = (double)compute_rgr(fw0, in.growth.total_biomass, (real_t)(cfg.days-1));
        int spd = (int)(24.0/cfg.dt_hours);
        int base = (cfg.days-1)*spd;
        double dt = (double)cfg.dt_hours;

        /* Re-evaluate the algebraic layer at each sampled point of day 30 and
         * integrate the terms with the same rectangle rule the sampler uses. */
        double Ia=0,Isd=0,Iu=0,It=0,Im=0,Ie=0,Iv=0;
        for (int k=0; k<spd && base+k<r.filled_steps; ++k) {
            double zt = k*dt;
            Input s = in;
            s.carbohydrates.sucrose         = r.sucrose[base+k];
            s.carbohydrates.starch          = r.starch[base+k];
            s.carbohydrates.sucrose_vacuole = r.sucrose_vacuole[base+k];
            /* k_s (X_STARCH_PART) is NOT sampled into SimulationResult -- the
             * buffer that looked like it was is total biomass.  Leave the
             * partition coefficient at its end-of-run value rather than
             * feeding biomass in as a fraction. */
            /* starch_night_start: the starch level at the last dusk, which is
             * what the running simulation would have latched. */
            int dusk = base + (int)(PP[i]/dt);
            if (dusk >= r.filled_steps) dusk = r.filled_steps-1;
            simulation_compute_algebraic(&s, (real_t)zt, r.starch[dusk]);
            double sp = (double)s.carbohydrates.starch_partition_coeff;
            Ia  += (1.0-sp)*(double)s.photo.photosynthesis*dt;
            Isd += (double)s.carbohydrates.starch_degradation_rate*dt;
            Iu  += (double)s.carbohydrates.uptake_cost*dt;
            It  += (double)s.carbohydrates.transport_cost*dt;
            /* The pre-change tree has no maintenance_respiration field: the
             * term was computed inline as respiration_frequency * sucrose, in
             * every phase.  Compiling this tool against both trees is how the
             * before/after rows are made comparable, so select the expression
             * at compile time rather than forking the tool. */
#ifdef PLANTMODEL_PRE_RESPIRATION_FIX
            Im  += (double)s.carbohydrates.respiration_frequency
                   *(double)s.carbohydrates.sucrose*dt;
#else
            Im  += (double)(s.carbohydrates.respiration_frequency * s.carbohydrates.sucrose)*dt;
#endif
            Ie  += (double)s.carbohydrates.sucrose_loading_frequency
                   *(double)s.carbohydrates.night_efficiency_starch*dt;
            Iv  += ((double)s.carbohydrates.k_vacuole_in*(double)s.carbohydrates.sucrose
                   -(double)s.carbohydrates.k_vacuole_out*(double)s.carbohydrates.sucrose_vacuole)*dt;
        }
        int last = base+spd < r.filled_steps ? base+spd : r.filled_steps-1;
        double dpool = (double)r.sucrose[last]-(double)r.sucrose[base];
        double resid = Ia+Isd-Iu-It-Im-Ie-Iv-dpool;
        /* Carbon locked into new biomass over day 30, on the measured basis:
         * RGR * FW * carbon density.  Reported as the density the model would
         * need for the export flux to be the growth flux. */
        double bC = Ie;                       /* what the model actually removes */
        double cdens = rgr > 0 ? bC/rgr : NAN;
        printf("%3dh %10.2f %9.2f %8.2f %8.2f %9.2f %10.2f %8.2f | %9.3f %10.3f %10.2f\n",
               PP[i], Ia, Isd, Iu, It, Im, Ie, Iv, dpool, resid, bC);
        fprintf(fp, "%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.4f,%.8f\n",
                PP[i], Ia, Isd, Iu, It, Im, Ie, Iv, dpool, resid, bC, cdens, rgr);
        simulation_result_free(&r);
    }
    fclose(fp);
}

int main(int argc, char** argv)
{
    const char* mode = "ceiling";
    const char* outdir = "results/recalibration";
    int iters = 90;
    double pfix[NPAR] = {0.0045, 2.0, 0.15, 0.08, 1.98};
    int have_pfix = 0;

    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i],"--mode") && i+1<argc)        mode = argv[++i];
        else if (!strcmp(argv[i],"--out") && i+1<argc)    outdir = argv[++i];
        else if (!strcmp(argv[i],"--iters") && i+1<argc)  iters = atoi(argv[++i]);
        else if (!strcmp(argv[i],"--release-clamp"))      g_release_clamp = 1;
        else if (!strcmp(argv[i],"--params") && i+5<argc) {
            for (int j=0;j<NPAR;++j) pfix[j] = atof(argv[++i]);
            have_pfix = 1;
        }
    }
    (void)have_pfix;

    if (!strcmp(mode,"ceiling"))      mode_ceiling(outdir, iters);
    else if (!strcmp(mode,"heldout")) mode_heldout(outdir, iters);
    else if (!strcmp(mode,"diel"))    mode_diel(outdir, pfix);
    else if (!strcmp(mode,"budget"))  mode_budget(outdir, pfix);
    else if (!strcmp(mode,"predict")) mode_predict(pfix);
    else { fprintf(stderr,"unknown mode %s\n", mode); return 1; }
    return 0;
}
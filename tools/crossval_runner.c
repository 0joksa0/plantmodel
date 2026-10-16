/*
 * crossval_runner.c — out-of-sample validation of the photoperiod response
 *
 * Why this tool exists
 * --------------------
 * The manuscript reports a 0.24 % mean absolute percentage error between
 * simulated and published relative growth rates across the four photoperiods.
 * That number is not a measure of predictive accuracy.  In
 * tools/ablation_runner.c the growth coefficient lambda_sb is obtained by
 * bisection until the simulated RGR equals the published RGR, separately for
 * each photoperiod, and lambda_sdi is obtained by an outer bisection against
 * the observed end-of-day starch.  With four photoperiods that is eight fitted
 * parameters against eight observations: zero degrees of freedom.  The
 * reported "errors" are the convergence residuals of the bisection, which is
 * why re-running the calibration returns 0.00 % for every photoperiod.
 *
 * What is tested here
 * -------------------
 * A model that fits one parameter per condition cannot be said to predict the
 * response to that condition.  To obtain a falsifiable statement we impose a
 * photoperiod-INDEPENDENT parameterisation, fit it on a subset of the
 * photoperiods, freeze it, and evaluate the held-out photoperiods.  Two
 * parameterisations are examined:
 *
 *   H1 "constant"  lambda_sb = c0,             lambda_sdi = s0     (2 params)
 *   H2 "linear"    lambda_sb = c0 + c1 * P,    lambda_sdi = s0     (3 params)
 *
 * H2 encodes the hypothesis that the growth coefficient varies smoothly with
 * photoperiod; it is still a prediction for a held-out photoperiod because the
 * held-out observation never enters the fit.
 *
 * Splits: the three balanced 2-vs-2 partitions of {4, 6, 8, 12} h, plus the
 * four leave-one-out partitions.  For each split the training objective is the
 * relative sum of squares over the training observables, minimised by
 * Nelder-Mead; predictions on the held-out photoperiods are then reported as
 * signed relative error, MAPE and nRMSE.
 *
 * A "saturated" reference row reproduces the manuscript's per-photoperiod
 * calibration so that the in-sample and out-of-sample numbers appear side by
 * side.
 *
 * Honest caveat, reported with the results: lambda_c (feedback on
 * photosynthesis) and lambda_sni (night starch partition) remain at their
 * per-photoperiod Table A.1 values in every configuration, because they are
 * taken from the source publication rather than fitted here.  The
 * out-of-sample test therefore isolates the growth and starch-synthesis
 * coefficients, and is a test of those, not of the whole parameter set.
 *
 * Output
 * ------
 *   results/crossval/crossval_summary.csv     one row per (hypothesis, split, photoperiod)
 *   results/crossval/crossval_params.csv      fitted parameters per (hypothesis, split)
 */

#define SOLVER_REAL_AS double
#define _POSIX_C_SOURCE 200809L

#include "config/config.h"
#include "model/input.h"
#include "model/model.h"
#include "simulation/simulation.h"

#include <math.h>
#include <solver.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

/* ------------------------------------------------------------------ observations */

/* Published values, identical to those used by tools/ablation_runner.c.
 *
 * SOURCE: all four photoperiod points come from Sulpice et al. (2014),
 * "Arabidopsis Coordinates the Diurnal Regulation of Carbon Allocation and
 * Growth across a Wide Range of Photoperiods", Mol Plant 7(1):137-155,
 * doi:10.1093/mp/sst127 (BibTeX key Sulpice2014).  That study grew Col-0
 * across a range of photoperiods and reported, per photoperiod:
 *   - rgr_obs    <- relative growth rates modelled from rosette biomass;
 *   - starch_obs <- rosette starch assayed at dusk (end of the light period),
 *                   i.e. ZT = photoperiod_h, which is what simulate_one()
 *                   samples.
 * Both columns therefore share a single source and a single growth series,
 * which is what makes the 4/6/8/12 h set internally consistent.
 *
 * Earlier revisions of this file credited Annunziata et al. (2017) and
 * Mengin et al. (2017).  That attribution was wrong for BOTH columns and
 * must not be reinstated: Annunziata et al. (2017) is "Getting back to
 * nature: a reality check for experiments in controlled environments"
 * (J Exp Bot 68(16):4463-4477), a sunlight-vs-artificial-light comparison
 * run at a 12 h photoperiod ONLY, and Mengin et al. (2017)
 * (Plant Cell Environ 40(11):2608-2627) covers 6 h and 12 h ONLY.  Neither
 * paper reports a 4 h or an 8 h point, so neither can be the source of this
 * table.  Those two studies remain legitimate references for diel
 * sucrose/starch *shape* elsewhere in the project, but not for these values.
 */
typedef struct {
    int    photoperiod_h;
    double rgr_obs;         /* published RGR [1/d]  (compute_rgr() is called
                             * with a duration in DAYS, not hours)          */
    double starch_obs;      /* starch at dusk [µmol/gFW]        */
    double lambda_c;        /* Table A.1, not fitted            */
    double lambda_sni;      /* Table A.1, not fitted            */
} Obs;

static const Obs OBS[] = {
    {  4, 0.0680, 29.67, 0.82, 0.16 },
    {  6, 0.1135, 54.24, 0.79, 0.15 },
    {  8, 0.1708, 56.88, 0.67, 0.13 },
    { 12, 0.2600, 58.73, 0.62, 0.08 },
};
#define NPP 4

/* ------------------------------------------------------------------ simulation */

static SimulationConfig g_cfg;
static Input            g_base;

/* `diverged` carries the simulation's own divergence verdict: a non-finite
 * trajectory, or an exhausted solver step budget (the adaptive step collapsing
 * without the state going non-finite — see SIMULATION_MAX_SOLVER_STEPS in
 * src/simulation/simulation.c).  Either way the reported RGR belongs to a
 * truncated trajectory and must not enter an objective or a prediction. */
typedef struct { double rgr, eod_starch; int diverged; } SimOut;

/* lambda_sni and lambda_c are passed explicitly rather than read from OBS[].
 *
 * In the original design both were taken from the per-photoperiod Table A.1
 * row, which meant that even the nominally two-parameter "constant" hypothesis
 * silently carried eight photoperiod-specific numbers.  A hypothesis claiming a
 * photoperiod-independent parameterisation cannot rest on eight per-photoperiod
 * constants, so the caller now decides: pass the Table A.1 values to reproduce
 * the earlier behaviour, or pass a single shared fitted value to test the
 * stronger claim.  Passing a negative value selects the Table A.1 default. */
static SimOut simulate_one_full(int pp_index, double lambda_sb, double lambda_sdi,
                                double lambda_sni, double lambda_c)
{
    const Obs* o = &OBS[pp_index];
    Input in = g_base;
    in.photo.photo_mode                  = (int)PHOTO_MODE_FARQUHAR_JARVIS;
    in.photo.photoperiod                 = (real_t)o->photoperiod_h;
    in.photo.feedback_on_photosynthesis  = (real_t)((lambda_c   < 0.0) ? o->lambda_c   : lambda_c);
    in.carbohydrates.lambda_sni          = (real_t)((lambda_sni < 0.0) ? o->lambda_sni : lambda_sni);
    in.carbohydrates.lambda_sdi          = (real_t)lambda_sdi;
    in.growth.lambda_sb                  = (real_t)lambda_sb;

    real_t fw_start = in.growth.total_biomass;

    SimulationResult r;
    memset(&r, 0, sizeof(r));
    simulate_days(&g_cfg, &in, &r);

    SimOut out;
    out.rgr = (double)compute_rgr(fw_start, in.growth.total_biomass,
                                  (real_t)(g_cfg.days - 1));

    int steps_per_day = (int)(REAL(24.0) / g_cfg.dt_hours);
    int last_start    = r.filled_steps - steps_per_day;
    if (last_start < 0) last_start = 0;
    int eod = last_start + (int)((real_t)o->photoperiod_h / g_cfg.dt_hours);
    if (eod >= r.filled_steps) eod = r.filled_steps - 1;
    out.eod_starch = (r.filled_steps > 0) ? (double)r.starch[eod] : 0.0;
    out.diverged   = r.diverged;

    simulation_result_free(&r);
    return out;
}

/* Table A.1 defaults for lambda_sni and lambda_c. */
static SimOut simulate_one(int pp_index, double lambda_sb, double lambda_sdi)
{
    return simulate_one_full(pp_index, lambda_sb, lambda_sdi, -1.0, -1.0);
}

/* ------------------------------------------------------------------ hypotheses */

/* The ladder is ordered by what the hypothesis is allowed to know about the
 * photoperiod, and the free-parameter count is stated honestly: a parameter
 * frozen at a per-photoperiod literature value is still information about the
 * photoperiod, even when it is not fitted.
 *
 *   H1  shared lambda_sb,          shared lambda_sni, shared lambda_c   (4)
 *   H2  lambda_sb linear in P,     shared lambda_sni, shared lambda_c   (5)
 *   H3  shared lambda_sb,          per-photoperiod lambda_sni/lambda_c  (2 + 8 frozen)
 *   H4  lambda_sb linear in P,     per-photoperiod lambda_sni/lambda_c  (3 + 8 frozen)
 *   H5  saturated per-photoperiod everything                            (8)
 *
 * H3 and H4 are the previously published "constant" and "linear" cases and are
 * retained unchanged for comparability.  H1 against H3 is the measurement this
 * ladder exists to make: it isolates how much of the apparent predictive skill
 * came from the eight frozen per-photoperiod numbers. */
typedef enum {
    H_SHARED_CONST  = 0,
    H_SHARED_LINEAR = 1,
    H_CONST         = 2,
    H_LINEAR        = 3,
    H_SATURATED     = 4
} Hypothesis;

#define NHYP 4   /* cross-validated hypotheses; H_SATURATED is in-sample only */

static const char* H_NAME[5] = {
    "shared_constant", "shared_linear",
    "constant", "linear", "saturated_per_photoperiod"
};
/* free parameters actually fitted */
static const int   H_NPAR[5] = { 4, 5, 2, 3, 8 };
/* additional per-photoperiod values frozen at Table A.1 (not fitted, but still
 * photoperiod-specific information supplied to the model) */
static const int   H_NFROZEN_PP[5] = { 0, 0, 8, 8, 0 };
/* dimension of the optimiser's parameter vector */
static const int   H_NDIM[5] = { 4, 5, 2, 3, 0 };

/* Map hypothesis parameter vector -> (lambda_sb, lambda_sdi) for a photoperiod.
 * Parameters are handled in a transformed space so the optimiser is
 * unconstrained: lambda_sdi and the constant term are exponentiated. */
static void unpack_params(Hypothesis h, const double* p, int pp_index,
                          double* lambda_sb, double* lambda_sdi,
                          double* lambda_sni, double* lambda_c)
{
    double P = (double)OBS[pp_index].photoperiod_h;

    /* negative => simulate_one_full() substitutes the Table A.1 value */
    *lambda_sni = -1.0;
    *lambda_c   = -1.0;

    switch (h) {
    case H_CONST:                       /* p = [log lsb, log lsdi] */
        *lambda_sb  = exp(p[0]);
        *lambda_sdi = exp(p[1]);
        break;
    case H_LINEAR: {                    /* p = [log c0, c1, log lsdi] */
        double v = exp(p[0]) + p[1] * P;
        if (v < 1e-6) v = 1e-6;
        *lambda_sb  = v;
        *lambda_sdi = exp(p[2]);
        break;
    }
    case H_SHARED_CONST:                /* p = [log lsb, log lsdi, logit sni, logit c] */
        *lambda_sb  = exp(p[0]);
        *lambda_sdi = exp(p[1]);
        /* lambda_sni and lambda_c are fractions; a logistic map keeps the
         * optimiser unconstrained while holding them inside (0,1). */
        *lambda_sni = 1.0 / (1.0 + exp(-p[2]));
        *lambda_c   = 1.0 / (1.0 + exp(-p[3]));
        break;
    case H_SHARED_LINEAR: {             /* p = [log c0, c1, log lsdi, logit sni, logit c] */
        double v = exp(p[0]) + p[1] * P;
        if (v < 1e-6) v = 1e-6;
        *lambda_sb  = v;
        *lambda_sdi = exp(p[2]);
        *lambda_sni = 1.0 / (1.0 + exp(-p[3]));
        *lambda_c   = 1.0 / (1.0 + exp(-p[4]));
        break;
    }
    default:
        *lambda_sb = *lambda_sdi = 0.0;
        break;
    }
}

/* Upper admissible bound on lambda_sb.
 *
 * This is a feasibility bound, not a tuning knob, and it exists for a
 * numerical reason.  The integrator cost grows linearly with lambda_sb once
 * the biomass trajectory becomes stiff: measured wall time for one 30-day
 * simulation is 12 ms at the calibrated scale (~5e-3), 0.5 s at 1e3, 4.9 s at
 * 1e4, 49 s at 1e5 and 414 s at 1e6.  The run never diverges and RGR stays
 * finite (0.86 to 1.09 d-1), so the isfinite() rejection below never fires and
 * a single unlucky simplex excursion can stall a fit for hours.
 *
 * The bound is placed far above the feasible region: the calibrated values are
 * 3.8e-3 to 5.9e-3, and any lambda_sb of order unity already produces RGR
 * above 0.8 d-1, more than three times the largest observation (0.26 d-1 at
 * 12 h).  Rejecting before the simulation runs therefore removes the runtime
 * pathology without excluding any parameterisation that could fit the data. */
#define LAMBDA_SB_MAX 1.0

/* Upper admissible bound on a SHARED lambda_sdi, for the same reason and of the
 * same kind as LAMBDA_SB_MAX above: a numerical feasibility bound, not a tuning
 * knob.
 *
 * Above a measured threshold of lambda_sdi ~ 0.5163 the adaptive step collapses
 * at the 12 h photoperiod and the run exhausts the 10 M step budget without the
 * state ever becoming non-finite.  The simulation is then correctly classified
 * as divergent, but only after minutes of wall time, and a single simplex
 * excursion into that region stalls the whole fit -- which is exactly what
 * happened to five shared_linear fits before this bound was added.
 *
 * The bound applies to the cross-validated hypotheses, which impose ONE
 * lambda_sdi across all four photoperiods and therefore must remain integrable
 * at 12 h.  It does not apply to fit_saturated(), which fits lambda_sdi
 * separately per photoperiod inside its own bracket and legitimately reaches
 * 0.60 at 6 h, where no stall exists. */
#define LAMBDA_SDI_MAX 0.50

/* Admissibility of a parameter vector, checked over ALL four photoperiods
 * rather than only the training ones.
 *
 * Each hypothesis asserts a photoperiod-INDEPENDENT law valid across the whole
 * 4-12 h range, so a coefficient pair that implies an inadmissible lambda_sb at
 * any photoperiod is rejected even when it is well behaved on the training
 * subset.  This matters for the linear hypothesis fitted on a narrow lever arm
 * (e.g. training on 4 and 6 h only): the slope is poorly determined there, and
 * extrapolation to 12 h can imply a lambda_sb several orders of magnitude above
 * the feasible scale.  Without this check the optimiser is free to accept such
 * a vector, and the subsequent prediction step then spends minutes integrating
 * a stiff trajectory that is physically meaningless. */
static int params_admissible(Hypothesis h, const double* p)
{
    for (int i = 0; i < NPP; ++i) {
        double lsb, lsdi, lsni, lc;
        unpack_params(h, p, i, &lsb, &lsdi, &lsni, &lc);
        if (!(lsb > 0.0) || !(lsdi > 0.0) || lsdi > LAMBDA_SDI_MAX) return 0;
        if (lsb > LAMBDA_SB_MAX) return 0;
        if (lsni >= 0.0 && !(lsni > 0.0 && lsni < 1.0)) return 0;
        if (lc   >= 0.0 && !(lc   > 0.0 && lc   < 1.0)) return 0;
    }
    return 1;
}

/* Relative sum of squares over a training mask. */
static double objective(Hypothesis h, const double* p, const int* mask)
{
    double sse = 0.0;
    if (!params_admissible(h, p)) return 1e30;
    for (int i = 0; i < NPP; ++i) {
        if (!mask[i]) continue;
        double lsb, lsdi, lsni, lc;
        unpack_params(h, p, i, &lsb, &lsdi, &lsni, &lc);
        SimOut s = simulate_one_full(i, lsb, lsdi, lsni, lc);
        if (s.diverged) return 1e30;
        if (!isfinite(s.rgr) || !isfinite(s.eod_starch)) return 1e30;
        double e1 = (s.rgr        - OBS[i].rgr_obs)    / OBS[i].rgr_obs;
        double e2 = (s.eod_starch - OBS[i].starch_obs) / OBS[i].starch_obs;
        sse += e1 * e1 + e2 * e2;
    }
    return sse;
}

/* ------------------------------------------------------------------ Nelder-Mead */

/* Largest optimiser dimension in the ladder (H_SHARED_LINEAR uses 5). */
#define NM_MAXDIM 5

typedef struct { Hypothesis h; const int* mask; } NMCtx;

static double nm_f(const double* x, void* vctx)
{
    NMCtx* c = (NMCtx*)vctx;
    return objective(c->h, x, c->mask);
}

static void nelder_mead(int n, double* x0, double step,
                        double (*f)(const double*, void*), void* ctx,
                        int max_iter, double tol, double* out_x, double* out_f)
{
    int    np = n + 1;
    double simplex[NM_MAXDIM + 1][NM_MAXDIM];
    double fv[NM_MAXDIM + 1];

    for (int i = 0; i < np; ++i) {
        for (int k = 0; k < n; ++k) simplex[i][k] = x0[k];
        if (i > 0) simplex[i][i - 1] += step;
        fv[i] = f(simplex[i], ctx);
    }

    for (int iter = 0; iter < max_iter; ++iter) {
        /* order: find best, worst, second-worst */
        int b = 0, w = 0, sw = 0;
        for (int i = 1; i < np; ++i) { if (fv[i] < fv[b]) b = i; if (fv[i] > fv[w]) w = i; }
        sw = (w == 0) ? 1 : 0;
        for (int i = 0; i < np; ++i) if (i != w && fv[i] > fv[sw]) sw = i;

        if (fabs(fv[w] - fv[b]) <= tol * (fabs(fv[b]) + fabs(fv[w]) + 1e-12))
            break;

        double centroid[NM_MAXDIM] = {0};
        for (int i = 0; i < np; ++i) {
            if (i == w) continue;
            for (int k = 0; k < n; ++k) centroid[k] += simplex[i][k];
        }
        for (int k = 0; k < n; ++k) centroid[k] /= (double)n;

        double xr[NM_MAXDIM], fr;
        for (int k = 0; k < n; ++k) xr[k] = centroid[k] + 1.0 * (centroid[k] - simplex[w][k]);
        fr = f(xr, ctx);

        if (fr < fv[b]) {                                  /* expansion */
            double xe[NM_MAXDIM], fe;
            for (int k = 0; k < n; ++k) xe[k] = centroid[k] + 2.0 * (centroid[k] - simplex[w][k]);
            fe = f(xe, ctx);
            if (fe < fr) { for (int k = 0; k < n; ++k) simplex[w][k] = xe[k]; fv[w] = fe; }
            else          { for (int k = 0; k < n; ++k) simplex[w][k] = xr[k]; fv[w] = fr; }
        } else if (fr < fv[sw]) {                          /* reflection */
            for (int k = 0; k < n; ++k) simplex[w][k] = xr[k];
            fv[w] = fr;
        } else {                                           /* contraction */
            double xc[NM_MAXDIM], fc;
            for (int k = 0; k < n; ++k) xc[k] = centroid[k] + 0.5 * (simplex[w][k] - centroid[k]);
            fc = f(xc, ctx);
            if (fc < fv[w]) { for (int k = 0; k < n; ++k) simplex[w][k] = xc[k]; fv[w] = fc; }
            else {                                          /* shrink */
                for (int i = 0; i < np; ++i) {
                    if (i == b) continue;
                    for (int k = 0; k < n; ++k)
                        simplex[i][k] = simplex[b][k] + 0.5 * (simplex[i][k] - simplex[b][k]);
                    fv[i] = f(simplex[i], ctx);
                }
            }
        }
    }

    int b = 0;
    for (int i = 1; i < np; ++i) if (fv[i] < fv[b]) b = i;
    for (int k = 0; k < n; ++k) out_x[k] = simplex[b][k];
    *out_f = fv[b];
}

/* Per-photoperiod saturated fit: bisection on lambda_sb for RGR, nested in a
 * bisection on lambda_sdi for end-of-day starch.  Reproduces the manuscript
 * procedure so its in-sample residual can be reported alongside. */
static void fit_saturated(int i, double* out_lsb, double* out_lsdi)
{
    double lo_sdi = 0.02, hi_sdi = 0.60;

    for (int it = 0; it < 40; ++it) {
        double mid_sdi = 0.5 * (lo_sdi + hi_sdi);

        double lo_sb = 1e-4, hi_sb = 0.10, mid_sb = 0.0;
        for (int j = 0; j < 60; ++j) {
            mid_sb = 0.5 * (lo_sb + hi_sb);
            double r = simulate_one(i, mid_sb, mid_sdi).rgr;
            if (r < OBS[i].rgr_obs) lo_sb = mid_sb; else hi_sb = mid_sb;
            if (hi_sb - lo_sb < 1e-10) break;
        }
        double eod = simulate_one(i, mid_sb, mid_sdi).eod_starch;
        if (eod < OBS[i].starch_obs) lo_sdi = mid_sdi; else hi_sdi = mid_sdi;
        *out_lsb  = mid_sb;
        *out_lsdi = mid_sdi;
        if (hi_sdi - lo_sdi < 1e-8) break;
    }
}

/* ------------------------------------------------------------------ splits */

typedef struct {
    const char* name;
    int train[NPP];   /* 1 = training, 0 = held out */
} Split;

static const Split SPLITS[] = {
    { "train_8_12__test_4_6",  { 0, 0, 1, 1 } },
    { "train_4_6__test_8_12",  { 1, 1, 0, 0 } },
    { "train_4_12__test_6_8",  { 1, 0, 0, 1 } },
    { "loo_holdout_4",         { 0, 1, 1, 1 } },
    { "loo_holdout_6",         { 1, 0, 1, 1 } },
    { "loo_holdout_8",         { 1, 1, 0, 1 } },
    { "loo_holdout_12",        { 1, 1, 1, 0 } },
};
#define NSPLIT ((int)(sizeof(SPLITS) / sizeof(SPLITS[0])))

/* Deterministic restarts per fit.  Eight covers the 4x2 grid of lambda_sb and
 * lambda_sdi seeds; the shared-fraction and slope seeds vary with the same
 * index so no two starts coincide. */
#define NSTART 8

/* 4 saturated reference fits + NHYP hypotheses x NSPLIT splits. */
#define NTASK (4 + NHYP * NSPLIT)

/* ------------------------------------------------------------------ main */

static void ensure_dir(const char* p) { mkdir(p, 0775); }

int main(int argc, char** argv)
{
    const char* config_path = "config/default.config";
    const char* out_dir     = "results/crossval";
    int         days_override = 0;
    /* --task k selects a single independent unit of work so that the 18 fits can
     * be run as concurrent processes instead of one serial sweep.  Each fit is a
     * separate Nelder-Mead (or bisection) optimisation over its own training
     * mask, so they share no state and the partition is exact rather than
     * approximate.  Tasks 0..3 are the saturated per-photoperiod reference fits,
     * 4..10 are hypothesis "constant" over the 7 splits, 11..17 are hypothesis
     * "linear" over the same splits.  Each task writes a headerless shard; the
     * shards are concatenated afterwards and the MAPE aggregates recomputed from
     * the merged rows, which is identical to computing them in one process. */
    int         task = -1;

    for (int i = 1; i < argc; ++i) {
        if      (strcmp(argv[i], "--config") == 0 && i + 1 < argc) config_path = argv[++i];
        else if (strcmp(argv[i], "--out")    == 0 && i + 1 < argc) out_dir     = argv[++i];
        else if (strcmp(argv[i], "--days")   == 0 && i + 1 < argc) days_override = atoi(argv[++i]);
        else if (strcmp(argv[i], "--task")   == 0 && i + 1 < argc) task = atoi(argv[++i]);
        else { fprintf(stderr, "Usage: %s [--config p] [--out d] [--days N] [--task k]\n", argv[0]); return 1; }
    }
    if (task >= NTASK) { fprintf(stderr, "--task must be in 0..%d (or omitted to run all)\n", NTASK - 1); return 1; }

    /* Start from the built-in defaults, exactly as ablation_runner does, then
     * overlay the config file.  config_load_file() only assigns the keys it
     * finds, so loading into a zero-initialised Input leaves every unlisted
     * field at 0 — which silently changes the model rather than failing. */
    g_cfg  = simulation_default_config();
    g_base = generate_input();

    if (config_load_file(config_path, &g_cfg, &g_base) != 0) {
        fprintf(stderr, "Failed to load config: %s\n", config_path);
        return 1;
    }
    g_cfg.gui_enabled = 0;
    g_cfg.write_csv   = 0;
    if (days_override > 0) g_cfg.days = days_override;

    ensure_dir("results");
    ensure_dir(out_dir);

    char p_summary[512], p_params[512], p_multi[512];
    if (task < 0) {
        snprintf(p_summary, sizeof(p_summary), "%s/crossval_summary.csv", out_dir);
        snprintf(p_params,  sizeof(p_params),  "%s/crossval_params.csv",  out_dir);
        snprintf(p_multi,   sizeof(p_multi),   "%s/crossval_multistart.csv", out_dir);
    } else {
        snprintf(p_summary, sizeof(p_summary), "%s/shard_%02d_summary.csv", out_dir, task);
        snprintf(p_params,  sizeof(p_params),  "%s/shard_%02d_params.csv",  out_dir, task);
        snprintf(p_multi,   sizeof(p_multi),   "%s/shard_%02d_multistart.csv", out_dir, task);
    }
    FILE* fs = fopen(p_summary, "w");
    FILE* fp = fopen(p_params,  "w");
    FILE* fm = fopen(p_multi,   "w");
    if (!fs || !fp || !fm) { fprintf(stderr, "cannot open output files\n"); return 1; }

    if (task < 0) {
        fprintf(fs, "hypothesis,n_params,split,photoperiod_h,role,"
                    "lambda_sb,lambda_sdi,rgr_obs,rgr_sim,rgr_rel_err_pct,"
                    "starch_obs,starch_sim,starch_rel_err_pct\n");
        fprintf(fp, "hypothesis,n_fitted,n_frozen_pp,split,train_sse,c0,c1,"
                    "lambda_sdi,shared_lambda_sni,shared_lambda_c\n");
        fprintf(fm, "hypothesis,n_fitted,split,n_started,n_converged,"
                    "sse_best,sse_worst,lambda_sb4_min,lambda_sb4_max,"
                    "lambda_sb4_spread_pct\n");
    }

    printf("=== Out-of-sample validation of the photoperiod response ===\n");
    printf("  config: %s   days: %d\n", config_path, g_cfg.days);
    printf("  observations: 4 photoperiods x {RGR, end-of-day starch} = 8\n\n");

    /* ---------------- saturated reference (the manuscript's procedure) ---------------- */
    printf("-- saturated: one lambda_sb and one lambda_sdi per photoperiod (8 params, 8 observations) --\n");
    double sat_mape_rgr = 0.0, sat_mape_st = 0.0;
    for (int i = 0; i < NPP; ++i) {
        if (task >= 0 && task != i) continue;
        double lsb, lsdi;
        fit_saturated(i, &lsb, &lsdi);
        SimOut s = simulate_one(i, lsb, lsdi);
        double e1 = 100.0 * (s.rgr - OBS[i].rgr_obs) / OBS[i].rgr_obs;
        double e2 = 100.0 * (s.eod_starch - OBS[i].starch_obs) / OBS[i].starch_obs;
        sat_mape_rgr += fabs(e1) / NPP;
        sat_mape_st  += fabs(e2) / NPP;
        printf("   %2dh  lambda_sb=%.8f  lambda_sdi=%.6f   RGR err=%+8.4f %%   starch err=%+8.4f %%\n",
               OBS[i].photoperiod_h, lsb, lsdi, e1, e2);
        fprintf(fs, "%s,%d,%s,%d,in_sample,%.10g,%.10g,%.6g,%.10g,%.6e,%.6g,%.10g,%.6e\n",
                H_NAME[H_SATURATED], H_NPAR[H_SATURATED], "none_all_fitted",
                OBS[i].photoperiod_h, lsb, lsdi,
                OBS[i].rgr_obs, s.rgr, e1, OBS[i].starch_obs, s.eod_starch, e2);
    }
    printf("   in-sample MAPE:  RGR %.4f %%   starch %.4f %%   (zero degrees of freedom)\n\n",
           sat_mape_rgr, sat_mape_st);

    /* ---------------- cross-validated hypotheses ---------------- */
    for (int hi = 0; hi < NHYP; ++hi) {
        Hypothesis h = (Hypothesis)hi;
        int n   = H_NDIM[hi];
        int nfp = H_NPAR[hi];
        printf("-- %s: %d fitted parameter(s), %d frozen per-photoperiod value(s) --\n",
               H_NAME[hi], nfp, H_NFROZEN_PP[hi]);

        for (int si = 0; si < NSPLIT; ++si) {
            /* Task ids: 0..3 saturated reference fits, then NHYP*NSPLIT
             * cross-validated fits laid out hypothesis-major:
             * task = 4 + hi*NSPLIT + si. */
            if (task >= 0 && task != 4 + hi * NSPLIT + si) continue;
            const Split* sp = &SPLITS[si];
            NMCtx ctx = { h, sp->train };

            /* MULTI-START.  A single simplex cannot tell a global optimum from
             * the first basin it falls into, and this objective runs a stiff
             * ODE inside every evaluation, so a flat or multimodal surface is a
             * real possibility rather than a hypothetical one.  Each fit is
             * therefore repeated from NSTART deterministic starting points
             * spanning the plausible range of each coefficient; the best is
             * reported and the spread across starts is written out as the
             * optimiser's own contribution to the uncertainty. */
            double xb[NM_MAXDIM], fb = 1e300;
            double start_f[NSTART];
            double start_lsb4[NSTART];   /* lambda_sb at 4 h, a comparable scalar */
            int    n_started = 0, n_converged = 0;

            for (int st = 0; st < NSTART; ++st) {
                double x0[NM_MAXDIM], xt[NM_MAXDIM], ft;
                /* deterministic ladder: lambda_sb spanning 1e-3..1e-2,
                 * lambda_sdi spanning 0.02..0.40, shared fractions near the
                 * Table A.1 range */
                double sb_seed  = 1.0e-3 * pow(10.0, (double)(st % 4) / 3.0);      /* 1e-3 .. 1e-2 */
                double sdi_seed = 0.02  * pow(20.0, (double)((st / 4) % 2));       /* 0.02, 0.40  */
                double sni_seed = (st % 2) ? 0.08 : 0.15;
                double lc_seed  = (st % 3) ? 0.70 : 0.62;
                double c1_seed  = ((st / 2) % 2) ? 2.0e-4 : -1.0e-4;

                switch (h) {
                case H_CONST:
                    x0[0] = log(sb_seed); x0[1] = log(sdi_seed); break;
                case H_LINEAR:
                    x0[0] = log(sb_seed); x0[1] = c1_seed; x0[2] = log(sdi_seed); break;
                case H_SHARED_CONST:
                    x0[0] = log(sb_seed); x0[1] = log(sdi_seed);
                    x0[2] = log(sni_seed / (1.0 - sni_seed));
                    x0[3] = log(lc_seed  / (1.0 - lc_seed));  break;
                case H_SHARED_LINEAR:
                    x0[0] = log(sb_seed); x0[1] = c1_seed; x0[2] = log(sdi_seed);
                    x0[3] = log(sni_seed / (1.0 - sni_seed));
                    x0[4] = log(lc_seed  / (1.0 - lc_seed));  break;
                default: continue;
                }

                nelder_mead(n, x0, 0.25, nm_f, &ctx, 400, 1e-10, xt, &ft);
                n_started++;
                start_f[st] = ft;
                {
                    double a, b2, c, d2;
                    unpack_params(h, xt, 0, &a, &b2, &c, &d2);
                    start_lsb4[st] = a;
                }
                if (ft < 1e29) n_converged++;
                if (ft < fb) { fb = ft; for (int k = 0; k < n; ++k) xb[k] = xt[k]; }
            }

            /* spread across starts: how much of the answer is the optimiser? */
            double f_best = fb, f_worst = -1e300, lsb_lo = 1e300, lsb_hi = -1e300;
            for (int st = 0; st < n_started; ++st) {
                if (start_f[st] > 1e29) continue;
                if (start_f[st] > f_worst) f_worst = start_f[st];
                if (start_lsb4[st] < lsb_lo) lsb_lo = start_lsb4[st];
                if (start_lsb4[st] > lsb_hi) lsb_hi = start_lsb4[st];
            }
            double lsb_spread_pct = (lsb_lo > 0.0 && lsb_hi > lsb_lo)
                                  ? 100.0 * (lsb_hi - lsb_lo) / lsb_lo : 0.0;
            fprintf(fm, "%s,%d,%s,%d,%d,%.8e,%.8e,%.10g,%.10g,%.4f\n",
                    H_NAME[hi], nfp, sp->name, n_started, n_converged,
                    f_best, (f_worst > -1e299 ? f_worst : NAN),
                    (lsb_lo < 1e299 ? lsb_lo : NAN),
                    (lsb_hi > -1e299 ? lsb_hi : NAN), lsb_spread_pct);

            double c0 = exp(xb[0]);
            int    lin = (h == H_LINEAR || h == H_SHARED_LINEAR);
            double c1  = lin ? xb[1] : 0.0;
            double sdi, sh_sni, sh_lc;
            {
                double a, b2, c, d2;
                unpack_params(h, xb, 0, &a, &b2, &c, &d2);
                sdi = b2; sh_sni = c; sh_lc = d2;
            }
            fprintf(fp, "%s,%d,%d,%s,%.8e,%.10g,%.10g,%.10g,%.10g,%.10g\n",
                    H_NAME[hi], nfp, H_NFROZEN_PP[hi], sp->name, fb, c0, c1, sdi,
                    sh_sni, sh_lc);

            double mape_in = 0.0, mape_out = 0.0, sse_out = 0.0;
            int n_in = 0, n_out = 0;
            printf("   %-22s ", sp->name);

            for (int i = 0; i < NPP; ++i) {
                double lsb, lsdi, lsni, lc;
                unpack_params(h, xb, i, &lsb, &lsdi, &lsni, &lc);
                /* Safety net mirroring params_admissible(): never integrate a
                 * trajectory outside the feasible scale, whatever the optimiser
                 * returned.  The cost is linear in lambda_sb (414 s per run at
                 * 1e6), so an unguarded call here can stall the whole fit. */
                if (lsb > LAMBDA_SB_MAX) {
                    fprintf(stderr,
                            "  [%s/%s] lambda_sb = %.6g at %d h exceeds the "
                            "admissible bound %g; recorded as non-admissible\n",
                            H_NAME[hi], sp->name, lsb, OBS[i].photoperiod_h,
                            (double)LAMBDA_SB_MAX);
                    fprintf(fs, "%s,%d,%s,%d,%s,%.10g,%.10g,%.6g,nan,nan,%.6g,nan,nan\n",
                            H_NAME[hi], nfp, sp->name, OBS[i].photoperiod_h,
                            sp->train[i] ? "train" : "test",
                            lsb, lsdi, OBS[i].rgr_obs, OBS[i].starch_obs);
                    continue;
                }
                SimOut s = simulate_one_full(i, lsb, lsdi, lsni, lc);
                if (s.diverged) {
                    fprintf(stderr,
                            "  [%s/%s] trajectory diverged at %d h "
                            "(lambda_sb = %.6g, lambda_sdi = %.6g); recorded as "
                            "non-admissible rather than scored\n",
                            H_NAME[hi], sp->name, OBS[i].photoperiod_h, lsb, lsdi);
                    fprintf(fs, "%s,%d,%s,%d,%s,%.10g,%.10g,%.6g,nan,nan,%.6g,nan,nan\n",
                            H_NAME[hi], nfp, sp->name, OBS[i].photoperiod_h,
                            sp->train[i] ? "train" : "test",
                            lsb, lsdi, OBS[i].rgr_obs, OBS[i].starch_obs);
                    continue;
                }
                double e1 = 100.0 * (s.rgr - OBS[i].rgr_obs) / OBS[i].rgr_obs;
                double e2 = 100.0 * (s.eod_starch - OBS[i].starch_obs) / OBS[i].starch_obs;

                fprintf(fs, "%s,%d,%s,%d,%s,%.10g,%.10g,%.6g,%.10g,%.6e,%.6g,%.10g,%.6e\n",
                        H_NAME[hi], nfp, sp->name, OBS[i].photoperiod_h,
                        sp->train[i] ? "train" : "test",
                        lsb, lsdi, OBS[i].rgr_obs, s.rgr, e1,
                        OBS[i].starch_obs, s.eod_starch, e2);

                if (sp->train[i]) { mape_in += fabs(e1); n_in++; }
                else {
                    mape_out += fabs(e1); n_out++;
                    double d = (s.rgr - OBS[i].rgr_obs) / OBS[i].rgr_obs;
                    sse_out += d * d;
                }
            }
            printf("train MAPE(RGR)=%7.3f %%   TEST MAPE(RGR)=%8.3f %%   test nRMSE=%7.4f\n",
                   mape_in / n_in, mape_out / n_out, sqrt(sse_out / n_out));
        }
        printf("\n");
    }

    fclose(fs);
    fclose(fp);
    fclose(fm);
    printf("  -> %s\n  -> %s\n  -> %s\n", p_summary, p_params, p_multi);
    printf("\n  Note: in the 'constant' and 'linear' hypotheses lambda_c and\n"
           "  lambda_sni remain at per-photoperiod Table A.1 values (8 frozen\n"
           "  photoperiod-specific numbers).  The 'shared_*' hypotheses fit one\n"
           "  value of each for all photoperiods, so their photoperiod-independence\n"
           "  claim is complete.\n");
    return 0;
}

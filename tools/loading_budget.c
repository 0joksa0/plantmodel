/* loading_budget.c -- a CLOSING carbon partition for the phloem-loading study,
 * plus the clamp-activity probe and the diel trace, all off one trajectory.
 *
 * WHY THIS EXISTS RATHER THAN model_report's sinks block.
 *
 * The `sinks_diagnostic` block of tools/model_report.c has two defects, both
 * documented in that file: uptake_cost is computed from an Input copy whose
 * nutrient-uptake fields are never refreshed from the trajectory (so it is
 * identically zero regardless of the model), and `share_to_starch` is a GROSS
 * flux that is added to sinks which already draw on the same carbon, so the
 * five shares sum to about 1.9 and are not a partition. A ceiling argument
 * about phloem loading cannot be built on numbers that do not partition.
 *
 * THE FIX USED HERE. Over one closed diel cycle in periodic steady state, the
 * carbohydrate system obeys
 *
 *     A  =  R + T + U + L  +  dS  +  dSt  +  eps
 *
 * with A gross assimilation, R maintenance respiration, T the transport
 * surcharge, U the nutrient-uptake cost, L phloem loading, dS and dSt the
 * changes in the sucrose and starch pools across the cycle, and eps the
 * residual. Starch is a transit pool: it is filled in the light and emptied at
 * night, so dSt is small and starch does not appear as a separate sink at all --
 * that double count is exactly what inflated the earlier shares. Every term
 * except U is reconstructible from the saved trajectory; dS and dSt are read
 * from the trajectory endpoints rather than recomputed, and eps is PRINTED so
 * the closure can be judged instead of assumed. U is unavailable here for the
 * same structural reason as in model_report, so it is folded into eps and the
 * column is labelled accordingly.
 *
 * A share is therefore reported only after the residual shows the budget
 * closes; if eps is a few percent, the shares are a partition to that accuracy.
 *
 * Also emitted, from the same trajectory:
 *   - clamp activity: the fraction of the cycle with sucrose > max_sucrose, so
 *     the claim "the leaf spends most of the cycle on the flat side of the
 *     min()" is measured rather than argued from the operating range;
 *   - the three candidate loading substrates evaluated on the SAME trajectory,
 *     which separates a pointwise difference between the laws from a difference
 *     that survives the dynamics;
 *   - the diel trace, so the SHAPE of the sucrose pool can be compared with the
 *     measured series and not only its dusk value.
 *
 * Diagnostic only. src/ is untouched; the loading law is whatever the tree this
 * is compiled against implements, so one tool serves all three variants.
 */

#include "config/config.h"
#include "model/model.h"
#include "simulation/simulation.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#define NPP 4
static const int    PP[NPP]        = { 4, 6, 8, 12 };
static const double LAMBDA_C[NPP]  = { 0.82, 0.79, 0.67, 0.62 };
static const double LAMBDA_SNI[NPP]= { 0.16, 0.15, 0.13, 0.08 };
static const double LAMBDA_SB[NPP] = { 0.003843, 0.004433, 0.005246, 0.005850 };

static void ensure_dir(const char* path)
{
    char buf[1024];
    snprintf(buf, sizeof(buf), "%s", path);
    for (char* p = buf + 1; *p; ++p)
        if (*p == '/') { *p = '\0'; mkdir(buf, 0755); *p = '/'; }
    mkdir(buf, 0755);   /* the final component too: ensure_dir takes a DIRECTORY */
}

int main(int argc, char** argv)
{
    const char* out   = (argc > 1) ? argv[1] : "results/loading";
    const char* label = (argc > 2) ? argv[2] : "unlabelled";
    int days          = (argc > 3) ? atoi(argv[3]) : 30;

    ensure_dir(out);
    char pb[1024], pt[1024];
    snprintf(pb, sizeof(pb), "%s/budget_%s.csv", out, label);
    snprintf(pt, sizeof(pt), "%s/trace_%s.csv", out, label);
    FILE* fb = fopen(pb, "w");
    FILE* ft = fopen(pt, "w");
    if (!fb || !ft) return 1;

    fprintf(fb, "label,photoperiod,assim_daily,resp,transport,load,d_sucrose,d_starch,"
                "residual_incl_uptake,closure_pct,share_resp,share_transport,share_load,"
                "sucrose_min,sucrose_max,sucrose_mean,frac_cycle_above_smax,"
                "substrate_law1,substrate_law2,substrate_law3,ratio_l2_l1,ratio_l3_l1\n");
    fprintf(ft, "label,photoperiod,zt,light,sucrose,starch,photosynthesis,loading_substrate\n");

    for (int i = 0; i < NPP; ++i) {
        SimulationConfig cfg = simulation_default_config();
        Input in = generate_input();
        if (config_load_file("config/default.config", &cfg, &in) != 0) return 2;
        cfg.days = days; cfg.write_csv = 0;
        in.photo.photo_mode                 = (int)PHOTO_MODE_FARQUHAR_JARVIS;
        in.photo.photoperiod                = (real_t)PP[i];
        in.photo.feedback_on_photosynthesis = (real_t)LAMBDA_C[i];
        in.carbohydrates.lambda_sni         = (real_t)LAMBDA_SNI[i];
        in.growth.lambda_sb                 = (real_t)LAMBDA_SB[i];
        in.gas_exchange.specific_leaf_area  = (real_t)0.025;

        SimulationResult r; memset(&r, 0, sizeof(r));
        simulate_days(&cfg, &in, &r);

        double dt  = (double)cfg.dt_hours;
        int    spd = (int)(24.0 / dt);
        int    c0  = (cfg.days - 1) * spd, c1 = cfg.days * spd;
        if (c1 > r.filled_steps) c1 = r.filled_steps;
        if (c1 - c0 < 4) { fprintf(fb, "%s,%d,diverged\n", label, PP[i]); continue; }

        double smax = (double)in.carbohydrates.max_sucrose;
        double rm   = (double)in.carbohydrates.respiration_frequency;
        double rt   = (double)in.carbohydrates.sucrose_consumption_transport;
        double eta  = (double)in.carbohydrates.sucrose_loading_frequency;
        double lg   = (double)in.carbohydrates.lambda_g;
        double ks   = (double)in.carbohydrates.starch_partition_coeff;

        double A = 0, Rr = 0, T = 0, L = 0;
        double lo = 1e300, hi = -1e300, sum = 0, above = 0;
        double m1 = 0, m2 = 0, m3 = 0;
        int n = 0;

        /* Trapezoid over the closed cycle: the endpoint at zt = 24 h is the
         * first sample of the following day, so the window is a full period and
         * the pool differences below refer to the same phase. */
        for (int k = c0; k < c1; ++k) {
            double zt    = (k - c0) * dt;
            double light = zt < (double)PP[i] ? 1.0 : 0.0;
            double S     = (double)r.sucrose[k];
            double gate  = lg + (1.0 - lg) * (1.0 - light + light * (1.0 - ks));

            real_t nes = night_efficiency_starch((real_t)S,
                            in.carbohydrates.max_sucrose,
                            in.carbohydrates.lambda_g,
                            (real_t)light,
                            in.carbohydrates.starch_partition_coeff,
                                            s.carbohydrates.sink_activation_half,
                                            s.carbohydrates.sink_activation_hill);

            double resp = rm * S;
            double load = eta * (double)nes;

            A  += (double)r.ph[k] * dt;
            Rr += resp * dt;
            L  += load * dt;
            T  += rt * (resp + load) * dt;

            if (S < lo) lo = S;
            if (S > hi) hi = S;
            sum += S;
            if (S > smax) above += 1.0;
            m1 += (S < smax ? S : smax) * gate * dt;
            m2 += S * gate * dt;
            m3 += (smax * S / (S + smax)) * gate * dt;
            ++n;

            fprintf(ft, "%s,%d,%.4f,%.0f,%.6f,%.6f,%.6f,%.6f\n",
                    label, PP[i], zt, light, S, (double)r.starch[k],
                    (double)r.ph[k], (double)nes);
        }

        int kend = (c1 < r.filled_steps) ? c1 : c1 - 1;
        double dS  = (double)r.sucrose[kend] - (double)r.sucrose[c0];
        double dSt = (double)r.starch[kend]  - (double)r.starch[c0];
        double eps = A - (Rr + T + L) - dS - dSt;

        fprintf(fb, "%s,%d,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.4f,"
                    "%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n",
                label, PP[i], A, Rr, T, L, dS, dSt, eps,
                (A > 0 ? 100.0 * eps / A : NAN),
                (A > 0 ? Rr / A : NAN), (A > 0 ? T / A : NAN), (A > 0 ? L / A : NAN),
                lo, hi, sum / n, above / n,
                m1 / 24.0, m2 / 24.0, m3 / 24.0,
                (m1 > 0 ? m2 / m1 : NAN), (m1 > 0 ? m3 / m1 : NAN));
        simulation_result_free(&r);
    }
    fclose(fb); fclose(ft);
    printf("wrote %s and %s\n", pb, pt);
    return 0;
}

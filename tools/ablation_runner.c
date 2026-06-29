#define SOLVER_REAL_AS double
#include "config/config.h"
#include "model/input.h"
#include "model/model.h"
#include "simulation/simulation.h"

#include <math.h>
#include <solver.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ types */

typedef struct {
    int    photoperiod_h;
    real_t lit_rgr;
    real_t lambda_c;               /* feedback_on_photosynthesis (Table A.1) */
    real_t lambda_sni;             /* night starch partition (Table A.1) */
    /*
     * Per-model canonical lambda_sb values:
     *   basic          — Table A.1 (Annunziata et al.), the reference model
     *   farquhar        — calibrated for corrected Farquhar C3 (gc=gs/1.6, J(PAR))
     *   farquhar_jarvis — calibrated for corrected Farquhar+Jarvis
     */
    real_t lambda_sb_basic;
    real_t lambda_sb_farquhar;
    real_t lambda_sb_farquhar_jarvis;
} Scenario;

static const Scenario SCENARIOS[] = {
/*  PP   lit_rgr   lc     lsni   lsb_basic  lsb_farq   lsb_fj    */
    {  4, REAL(0.0680), REAL(0.82), REAL(0.16), REAL(0.00344), REAL(0.004264), REAL(0.003843) },
    {  6, REAL(0.1135), REAL(0.79), REAL(0.15), REAL(0.00413), REAL(0.004852), REAL(0.004433) },
    {  8, REAL(0.1708), REAL(0.67), REAL(0.13), REAL(0.00515), REAL(0.005648), REAL(0.005246) },
    { 12, REAL(0.2600), REAL(0.62), REAL(0.08), REAL(0.00578), REAL(0.006112), REAL(0.005850) },
};

typedef struct {
    PhotosynthesisMode mode;
    const char*        name;
} ModeEntry;

static const ModeEntry MODES[] = {
    { PHOTO_MODE_BASIC,           "basic"           },
    { PHOTO_MODE_FARQUHAR,        "farquhar"        },
    { PHOTO_MODE_FARQUHAR_JARVIS, "farquhar_jarvis" },
};

typedef struct {
    real_t model_rgr;
    real_t rel_err_pct;
    real_t max_ph;
    real_t day_avg_ph;
    real_t night_avg_ph;
    real_t min_sucrose;
    real_t min_starch;
    real_t final_total_biomass;
    real_t lambda_sb_used;
    int    nonfinite_samples;
} RunMetrics;

/* ------------------------------------------------------------------ helpers */

static real_t scenario_lambda_sb(const Scenario* s, PhotosynthesisMode mode)
{
    switch (mode) {
    case PHOTO_MODE_FARQUHAR:        return s->lambda_sb_farquhar;
    case PHOTO_MODE_FARQUHAR_JARVIS: return s->lambda_sb_farquhar_jarvis;
    case PHOTO_MODE_BASIC:
    default:                         return s->lambda_sb_basic;
    }
}

static void apply_scenario(Input* in, const Input* base, const Scenario* s, PhotosynthesisMode mode)
{
    *in = *base;
    in->photo.photo_mode                 = (int)mode;
    in->photo.photoperiod                = (real_t)s->photoperiod_h;
    in->photo.feedback_on_photosynthesis = s->lambda_c;
    in->carbohydrates.lambda_sni         = s->lambda_sni;
    in->growth.lambda_sb                 = scenario_lambda_sb(s, mode);
}

static real_t eval_rgr(const SimulationConfig* cfg, Input* in, int days)
{
    real_t fw_start = in->growth.total_biomass;
    SimulationResult result;
    memset(&result, 0, sizeof(result));
    simulate_days(cfg, in, &result);
    simulation_result_free(&result);
    return compute_rgr(fw_start, in->growth.total_biomass, (real_t)(days - 1));
}

/* ------------------------------------------------------------------ calibration */

/*
 * Bisection on lambda_sb to hit target RGR.
 * lambda_c and lambda_sni stay at Table A.1 values — they govern day/night
 * dynamics and are fixed per the original paper.
 * pmax stays fixed (not a per-photoperiod parameter).
 */
static real_t calibrate_lambda_sb(
    const SimulationConfig* cfg,
    const Input*            base_input,
    const Scenario*         scenario,
    PhotosynthesisMode      mode,
    real_t                  target_rgr,
    int                     days)
{
    real_t lo = REAL(0.0001);
    real_t hi = REAL(0.10);

    Input in;

    apply_scenario(&in, base_input, scenario, mode);
    in.growth.lambda_sb = lo;
    real_t rgr_lo = eval_rgr(cfg, &in, days);

    apply_scenario(&in, base_input, scenario, mode);
    in.growth.lambda_sb = hi;
    real_t rgr_hi = eval_rgr(cfg, &in, days);

    if (target_rgr <= rgr_lo) return lo;
    if (target_rgr >= rgr_hi) return hi;

    for (int i = 0; i < 50; i++) {
        real_t mid = (lo + hi) * REAL(0.5);
        apply_scenario(&in, base_input, scenario, mode);
        in.growth.lambda_sb = mid;
        real_t rgr_mid = eval_rgr(cfg, &in, days);

        if (rgr_mid < target_rgr) lo = mid;
        else                       hi = mid;

        if ((hi - lo) < REAL(1e-8))
            break;
    }

    return (lo + hi) * REAL(0.5);
}

/* ------------------------------------------------------------------ run */

static void write_timeseries(const SimulationResult*, const char*, int, real_t);

static void run_scenario(
    const SimulationConfig* cfg,
    Input*                  input,
    real_t                  lit_rgr,
    int                     days,
    const char*             ts_mode_name,
    RunMetrics*             out)
{
    SimulationResult result;
    memset(&result, 0, sizeof(result));
    simulate_days(cfg, input, &result);

    real_t fw_start      = REAL(0.04);
    out->model_rgr       = compute_rgr(fw_start, input->growth.total_biomass, (real_t)(days - 1));
    out->rel_err_pct     = lit_rgr > REAL(0.0)
                           ? RABS(out->model_rgr - lit_rgr) / lit_rgr * REAL(100.0)
                           : REAL(0.0);
    out->lambda_sb_used      = input->growth.lambda_sb;
    out->max_ph              = REAL(0.0);
    out->day_avg_ph          = REAL(0.0);
    out->night_avg_ph        = REAL(0.0);
    out->min_sucrose         = result.filled_steps > 0 ? result.sucrose[0] : REAL(0.0);
    out->min_starch          = result.filled_steps > 0 ? result.starch[0]  : REAL(0.0);
    out->nonfinite_samples   = 0;
    out->final_total_biomass = input->growth.total_biomass;

    real_t day_sum = REAL(0.0), night_sum = REAL(0.0);
    int    day_cnt = 0,         night_cnt = 0;

    for (int k = 0; k < result.filled_steps; ++k) {
        real_t ph      = result.ph[k];
        real_t sucrose = result.sucrose[k];
        real_t starch  = result.starch[k];

        if (!isfinite(ph) || !isfinite(sucrose) || !isfinite(starch)) {
            out->nonfinite_samples++;
            continue;
        }
        if (ph > out->max_ph)           out->max_ph      = ph;
        if (sucrose < out->min_sucrose) out->min_sucrose = sucrose;
        if (starch  < out->min_starch)  out->min_starch  = starch;

        real_t h = fmod((real_t)k * cfg->dt_hours, REAL(24.0));
        if (h < input->photo.photoperiod) { day_sum   += ph; day_cnt++;   }
        else                              { night_sum += ph; night_cnt++; }
    }

    if (day_cnt   > 0) out->day_avg_ph   = day_sum   / (real_t)day_cnt;
    if (night_cnt > 0) out->night_avg_ph = night_sum / (real_t)night_cnt;

    if (ts_mode_name)
        write_timeseries(&result, ts_mode_name, (int)input->photo.photoperiod, cfg->dt_hours);

    simulation_result_free(&result);
}

/* ------------------------------------------------------------------ time series writer */

static void write_timeseries(
    const SimulationResult* result,
    const char*             mode_name,
    int                     photoperiod_h,
    real_t                  dt_hours)
{
    char path[512];
    snprintf(path, sizeof(path),
             "results/ablation/timeseries/%s_%dh.csv",
             mode_name, photoperiod_h);

    FILE* f = fopen(path, "w");
    if (!f) {
        fprintf(stderr, "Warning: could not write timeseries %s\n", path);
        return;
    }

    fprintf(f, "zt_h,sucrose,starch,photosynthesis\n");

    /* Export last full 24h day so dynamics are at steady state */
    int steps_per_day = (int)(REAL(24.0) / dt_hours);
    int last_day_start = result->filled_steps - steps_per_day;
    if (last_day_start < 0) last_day_start = 0;

    for (int k = last_day_start; k < result->filled_steps; ++k) {
        real_t zt = (real_t)(k - last_day_start) * dt_hours;
        fprintf(f, "%.2f,%.6f,%.6f,%.6f\n",
                (double)zt,
                (double)result->sucrose[k],
                (double)result->starch[k],
                (double)result->ph[k]);
    }

    fclose(f);
}

/* ------------------------------------------------------------------ config writer */

static void write_ablation_config(
    const char*        mode_name,
    const Scenario*    scenario,
    real_t             lambda_sb_calibrated,
    int                days,
    const char*        base_config_path)
{
    char path[512];
    snprintf(path, sizeof(path),
             "config/ablation/%s/%dh.config",
             mode_name, scenario->photoperiod_h);

    FILE* f = fopen(path, "w");
    if (!f) {
        fprintf(stderr, "Warning: could not write %s\n", path);
        return;
    }

    fprintf(f,
        "# Ablation config: %s model, photoperiod=%dh\n"
        "# Base config: %s\n"
        "#\n"
        "# Per-photoperiod parameters (Table A.1, Annunziata et al.):\n"
        "#   lambda_c   (feedback_on_photosynthesis) = %.4f  [fixed, not calibrated]\n"
        "#   lambda_sni (night starch partition)     = %.4f  [fixed, not calibrated]\n"
        "#   lambda_sb  (growth rate)                = %.8f  [calibrated to match lit RGR=%.4f]\n"
        "#\n"
        "# Reproduce: ./build/ablation_runner --days %d --config %s\n"
        "\n"
        "simulation.days = %d\n"
        "simulation.gui_enabled = false\n"
        "simulation.write_csv = false\n"
        "\n"
        "input.photo.photoperiod = %d\n"
        "input.photo.feedback_on_photosynthesis = %.6f\n"
        "input.carbohydrates.lambda_sni = %.6f\n"
        "input.growth.lambda_sb = %.8f\n",
        mode_name, scenario->photoperiod_h,
        base_config_path,
        (double)scenario->lambda_c,
        (double)scenario->lambda_sni,
        (double)lambda_sb_calibrated,
        (double)scenario->lit_rgr,
        days, base_config_path,
        days,
        scenario->photoperiod_h,
        (double)scenario->lambda_c,
        (double)scenario->lambda_sni,
        (double)lambda_sb_calibrated);

    fclose(f);
}

/* ------------------------------------------------------------------ main */

static void usage(const char* argv0)
{
    printf("Usage: %s [--config path] [--days N] [--out csv] [--no-calibrate]\n", argv0);
    printf("  --no-calibrate   Use Table A.1 lambda_sb values without calibration\n");
}

int main(int argc, char** argv)
{
    const char* config_path  = "config/default.config";
    const char* out_path     = "results/ablation/ablation_results.csv";
    int         days         = 30;
    int         do_calibrate = 1;

    for (int i = 1; i < argc; ++i) {
        if      (strcmp(argv[i], "--config")       == 0 && i + 1 < argc) config_path = argv[++i];
        else if (strcmp(argv[i], "--days")         == 0 && i + 1 < argc) days = atoi(argv[++i]);
        else if (strcmp(argv[i], "--out")          == 0 && i + 1 < argc) out_path = argv[++i];
        else if (strcmp(argv[i], "--no-calibrate") == 0) do_calibrate = 0;
        else if (strcmp(argv[i], "--help")         == 0 || strcmp(argv[i], "-h") == 0) {
            usage(argv[0]);
            return 0;
        } else {
            fprintf(stderr, "Unknown argument: %s\n", argv[i]);
            usage(argv[0]);
            return 1;
        }
    }

    if (days < 2) { fprintf(stderr, "--days must be >= 2\n"); return 1; }

    SimulationConfig base_cfg = simulation_default_config();
    Input base_input = generate_input();
    if (config_load_file(config_path, &base_cfg, &base_input) != 0) {
        fprintf(stderr, "Failed to load config: %s\n", config_path);
        return 1;
    }
    base_cfg.days        = days;
    base_cfg.gui_enabled = 0;
    base_cfg.write_csv   = 0;

    FILE* out = fopen(out_path, "w");
    if (!out) { fprintf(stderr, "Failed to open output: %s\n", out_path); return 1; }

    fprintf(out,
        "photo_mode,photoperiod_h,days,lit_rgr,model_rgr,rel_err_pct,"
        "lambda_c,lambda_sni,lambda_sb,"
        "max_ph,day_avg_ph,night_avg_ph,"
        "min_sucrose,min_starch,final_total_biomass,nonfinite_samples\n");

    int n_modes    = (int)(sizeof(MODES)     / sizeof(MODES[0]));
    int n_scenario = (int)(sizeof(SCENARIOS) / sizeof(SCENARIOS[0]));

    for (int m = 0; m < n_modes; ++m) {
        printf("\n=== Mode: %s ===\n", MODES[m].name);
        if (do_calibrate)
            printf("  Calibrating lambda_sb per photoperiod...\n");

        printf("  %-5s  %-8s  %-10s  %-9s  %-6s  %-6s  %-10s  %-10s\n",
               "PP_h", "lit_rgr", "model_rgr", "err%", "lc", "lsni", "lambda_sb", "final_BM");

        for (int s = 0; s < n_scenario; ++s) {
            real_t lsb = scenario_lambda_sb(&SCENARIOS[s], MODES[m].mode);

            if (do_calibrate) {
                lsb = calibrate_lambda_sb(
                    &base_cfg, &base_input, &SCENARIOS[s],
                    MODES[m].mode, SCENARIOS[s].lit_rgr, days);
            }

            Input in;
            apply_scenario(&in, &base_input, &SCENARIOS[s], MODES[m].mode);
            in.growth.lambda_sb = lsb;

            RunMetrics metrics;
            memset(&metrics, 0, sizeof(metrics));
            run_scenario(&base_cfg, &in, SCENARIOS[s].lit_rgr, days, MODES[m].name, &metrics);

            if (do_calibrate)
                write_ablation_config(MODES[m].name, &SCENARIOS[s], lsb, days, config_path);

            fprintf(out,
                "%s,%d,%d,%.6f,%.6f,%.3f,"
                "%.4f,%.4f,%.8f,"
                "%.6f,%.6f,%.6f,"
                "%.6f,%.6f,%.6f,%d\n",
                MODES[m].name,
                SCENARIOS[s].photoperiod_h, days,
                (double)SCENARIOS[s].lit_rgr,
                (double)metrics.model_rgr,
                (double)metrics.rel_err_pct,
                (double)SCENARIOS[s].lambda_c,
                (double)SCENARIOS[s].lambda_sni,
                (double)metrics.lambda_sb_used,
                (double)metrics.max_ph,
                (double)metrics.day_avg_ph,
                (double)metrics.night_avg_ph,
                (double)metrics.min_sucrose,
                (double)metrics.min_starch,
                (double)metrics.final_total_biomass,
                metrics.nonfinite_samples);

            printf("  %-5d  %-8.4f  %-10.4f  %-9.2f  %-6.2f  %-6.2f  %-10.6f  %-10.4f\n",
                   SCENARIOS[s].photoperiod_h,
                   (double)SCENARIOS[s].lit_rgr,
                   (double)metrics.model_rgr,
                   (double)metrics.rel_err_pct,
                   (double)SCENARIOS[s].lambda_c,
                   (double)SCENARIOS[s].lambda_sni,
                   (double)metrics.lambda_sb_used,
                   (double)metrics.final_total_biomass);
        }
    }

    fclose(out);

    if (do_calibrate)
        printf("\nCalibrated configs written to config/ablation/{basic,farquhar,farquhar_jarvis}/\n");
    printf("Ablation results written to %s\n", out_path);
    return 0;
}

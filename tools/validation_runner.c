#define SOLVER_REAL_AS double
#include "config/config.h"
#include "model/input.h"
#include "simulation/simulation.h"

#include <math.h>
#include <solver.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int photoperiod_h;
    real_t lit_rgr;
    real_t lambda_sni;
    real_t lambda_sb;
    real_t feedback;
} Scenario;

typedef struct {
    int photoperiod_h;
    real_t zt;
    real_t starch;
    real_t sucrose;
} ReferencePoint;

typedef struct {
    real_t model_rgr;
    real_t abs_err;
    real_t rel_err;
    real_t max_ph;
    real_t day_avg_ph;
    real_t night_avg_ph;
    real_t min_sucrose;
    real_t min_starch;
    real_t min_total_biomass;
    real_t max_biomass_drop;
    int negative_biomass_steps;
    int nonfinite_samples;
    real_t final_total_biomass;
    int ref_points;
    real_t sucrose_scale;
    real_t sucrose_nrmse;
    real_t starch_scale;
    real_t starch_nrmse;
} RunMetrics;

typedef struct {
    int pass_rgr;
    int pass_day_night;
    int pass_state;
    int pass_growth;
    int pass_dt;
    int pass_external;
    int pass;
    char fail_reason[256];
} PassFlags;

static const Scenario SCENARIOS[] = {
    { 4, REAL(0.0680), REAL(0.16), REAL(0.00344), REAL(0.82) },
    { 6, REAL(0.1135), REAL(0.15), REAL(0.00413), REAL(0.79) },
    { 8, REAL(0.1708), REAL(0.13), REAL(0.00515), REAL(0.67) },
    { 12, REAL(0.2600), REAL(0.08), REAL(0.00578), REAL(0.62) },
};

#define REF_POINTS_CAP 2048

static void usage(const char* argv0)
{
    printf("Usage: %s [--config path] [--days N] [--out csv] [--threshold-rel X] [--threshold-night-day-ratio X] [--min-growth-rel X] [--max-negative-step-abs X] [--threshold-dt-rgr-pct X] [--threshold-dt-biomass-pct X] [--sugars-starch-data path] [--threshold-sucrose-nrmse X] [--threshold-starch-nrmse X] [--min-ref-points N]\\n", argv0);
}

static void append_fail_reason(char* buf, size_t buflen, const char* reason)
{
    size_t len = strlen(buf);
    if (len >= buflen - 1)
        return;
    if (len > 0) {
        strncat(buf, "|", buflen - len - 1);
        len = strlen(buf);
    }
    strncat(buf, reason, buflen - len - 1);
}

static real_t rel_pct(real_t a, real_t b)
{
    real_t denom = RABS(a) > REAL(1e-12) ? RABS(a) : REAL(1e-12);
    return (RABS(a - b) / denom) * REAL(100.0);
}

static int parse_real_strict(const char* token, real_t* out)
{
    if (!token || !*token)
        return 0;

    char* end = NULL;
    double v = strtod(token, &end);
    if (end == token)
        return 0;
    if (!isfinite(v))
        return 0;

    while (*end == ' ' || *end == '\t' || *end == '\r' || *end == '\n')
        ++end;
    if (*end != '\0')
        return 0;

    *out = (real_t)v;
    return 1;
}

static int parse_photoperiod_hours(const char* token, int* out)
{
    if (!token || !*token)
        return 0;

    char* end = NULL;
    long h = strtol(token, &end, 10);
    if (end == token)
        return 0;

    while (*end == ' ' || *end == '\t')
        ++end;
    if (*end == 'h' || *end == 'H')
        ++end;
    while (*end == ' ' || *end == '\t' || *end == '\r' || *end == '\n')
        ++end;
    if (*end != '\0')
        return 0;

    if (h <= 0 || h > 24)
        return 0;

    *out = (int)h;
    return 1;
}

static int load_reference_data(const char* path, ReferencePoint* out, int cap, int* out_count)
{
    FILE* f = fopen(path, "r");
    if (!f)
        return -1;

    char line[2048];
    int count = 0;
    int line_no = 0;

    while (fgets(line, sizeof(line), f)) {
        line_no++;
        if (line_no == 1)
            continue;

        char* cols[16];
        int col_count = 0;
        char* save = NULL;
        for (char* tok = strtok_r(line, ",", &save); tok && col_count < 16; tok = strtok_r(NULL, ",", &save)) {
            cols[col_count++] = tok;
        }

        if (col_count < 7)
            continue;

        int photoperiod_h = 0;
        real_t zt = REAL(0.0);
        real_t starch = REAL(0.0);
        real_t sucrose = REAL(0.0);

        if (!parse_photoperiod_hours(cols[0], &photoperiod_h))
            continue;
        if (!parse_real_strict(cols[1], &zt))
            continue;
        if (!parse_real_strict(cols[3], &starch))
            continue;
        if (!parse_real_strict(cols[6], &sucrose))
            continue;

        if (count >= cap)
            break;

        out[count].photoperiod_h = photoperiod_h;
        out[count].zt = zt;
        out[count].starch = starch;
        out[count].sucrose = sucrose;
        count++;
    }

    fclose(f);
    *out_count = count;
    return 0;
}

static void fit_scaled_nrmse(const real_t* model, const real_t* ref, int n, real_t* out_scale, real_t* out_nrmse)
{
    if (n <= 0) {
        *out_scale = REAL(1.0);
        *out_nrmse = REAL(1e9);
        return;
    }

    real_t dot_mr = REAL(0.0);
    real_t dot_mm = REAL(0.0);
    real_t ref_mean = REAL(0.0);

    for (int i = 0; i < n; ++i) {
        dot_mr += model[i] * ref[i];
        dot_mm += model[i] * model[i];
        ref_mean += ref[i];
    }

    ref_mean /= (real_t)n;
    real_t scale = dot_mm > REAL(1e-12) ? dot_mr / dot_mm : REAL(1.0);

    real_t mse = REAL(0.0);
    for (int i = 0; i < n; ++i) {
        real_t e = scale * model[i] - ref[i];
        mse += e * e;
    }
    mse /= (real_t)n;

    real_t rmse = RSQRT(mse);
    real_t denom = RABS(ref_mean) > REAL(1e-12) ? RABS(ref_mean) : REAL(1.0);

    *out_scale = scale;
    *out_nrmse = rmse / denom;
}

static void compute_external_fit_metrics(const SimulationConfig* cfg, const SimulationResult* result,
    int photoperiod_h, const ReferencePoint* ref_data, int ref_count, RunMetrics* out)
{
    out->ref_points = 0;
    out->sucrose_scale = REAL(1.0);
    out->sucrose_nrmse = REAL(-1.0);
    out->starch_scale = REAL(1.0);
    out->starch_nrmse = REAL(-1.0);

    if (!ref_data || ref_count <= 0 || result->filled_steps <= 0)
        return;

    real_t model_sucrose[256];
    real_t ref_sucrose[256];
    real_t model_starch[256];
    real_t ref_starch[256];
    int n = 0;

    real_t last_day_start = (real_t)(cfg->days - 1) * REAL(24.0);

    for (int i = 0; i < ref_count && n < 256; ++i) {
        if (ref_data[i].photoperiod_h != photoperiod_h)
            continue;
        if (ref_data[i].zt < REAL(0.0) || ref_data[i].zt > REAL(24.0))
            continue;

        real_t t_target = last_day_start + ref_data[i].zt;
        int idx = (int)llround((double)(t_target / cfg->dt_hours));
        if (idx < 0)
            idx = 0;
        if (idx >= result->filled_steps)
            idx = result->filled_steps - 1;

        real_t m_suc = result->sucrose[idx];
        real_t m_sta = result->starch[idx];
        if (!isfinite(m_suc) || !isfinite(m_sta))
            continue;

        model_sucrose[n] = m_suc;
        ref_sucrose[n] = ref_data[i].sucrose;
        model_starch[n] = m_sta;
        ref_starch[n] = ref_data[i].starch;
        n++;
    }

    out->ref_points = n;
    if (n < 3)
        return;

    fit_scaled_nrmse(model_sucrose, ref_sucrose, n, &out->sucrose_scale, &out->sucrose_nrmse);
    fit_scaled_nrmse(model_starch, ref_starch, n, &out->starch_scale, &out->starch_nrmse);
}

static void compute_run_metrics(const SimulationConfig* cfg, Input* run_input, real_t fw_start, real_t lit_rgr,
    real_t max_negative_step_abs, const ReferencePoint* ref_data, int ref_count, RunMetrics* out)
{
    SimulationResult result;
    memset(&result, 0, sizeof(result));
    simulate_days(cfg, run_input, &result);

    real_t rgr_days = (real_t)(cfg->days - 1);
    real_t model_rgr = compute_rgr(fw_start, run_input->growth.total_biomass, rgr_days);

    out->model_rgr = model_rgr;
    out->abs_err = RABS(model_rgr - lit_rgr);
    out->rel_err = (lit_rgr > REAL(0.0)) ? out->abs_err / lit_rgr : REAL(0.0);

    out->max_ph = REAL(0.0);
    out->day_avg_ph = REAL(0.0);
    out->night_avg_ph = REAL(0.0);
    out->min_sucrose = REAL(0.0);
    out->min_starch = REAL(0.0);
    out->min_total_biomass = REAL(0.0);
    out->max_biomass_drop = REAL(0.0);
    out->negative_biomass_steps = 0;
    out->nonfinite_samples = 0;
    out->final_total_biomass = run_input->growth.total_biomass;

    out->ref_points = 0;
    out->sucrose_scale = REAL(1.0);
    out->sucrose_nrmse = REAL(-1.0);
    out->starch_scale = REAL(1.0);
    out->starch_nrmse = REAL(-1.0);

    if (result.filled_steps > 0) {
        out->max_ph = result.ph[0];
        out->min_sucrose = result.sucrose[0];
        out->min_starch = result.starch[0];
        out->min_total_biomass = result.partition[0];
    }

    real_t day_sum = REAL(0.0);
    real_t night_sum = REAL(0.0);
    int day_count = 0;
    int night_count = 0;

    for (int k = 0; k < result.filled_steps; ++k) {
        real_t ph = result.ph[k];
        real_t sucrose = result.sucrose[k];
        real_t starch = result.starch[k];
        real_t biomass = result.partition[k];

        if (!isfinite(ph) || !isfinite(sucrose) || !isfinite(starch) || !isfinite(biomass)) {
            out->nonfinite_samples++;
            continue;
        }

        if (ph > out->max_ph)
            out->max_ph = ph;
        if (sucrose < out->min_sucrose)
            out->min_sucrose = sucrose;
        if (starch < out->min_starch)
            out->min_starch = starch;
        if (biomass < out->min_total_biomass)
            out->min_total_biomass = biomass;

        real_t hour = (real_t)k * cfg->dt_hours;
        real_t hour_of_day = fmod(hour, REAL(24.0));
        if (hour_of_day < run_input->photo.photoperiod) {
            day_sum += ph;
            day_count++;
        } else {
            night_sum += ph;
            night_count++;
        }

        if (k > 0) {
            real_t d_biomass = biomass - result.partition[k - 1];
            if (d_biomass < -max_negative_step_abs) {
                out->negative_biomass_steps++;
                if (-d_biomass > out->max_biomass_drop)
                    out->max_biomass_drop = -d_biomass;
            }
        }
    }

    if (day_count > 0)
        out->day_avg_ph = day_sum / (real_t)day_count;
    if (night_count > 0)
        out->night_avg_ph = night_sum / (real_t)night_count;

    compute_external_fit_metrics(cfg, &result, (int)run_input->photo.photoperiod, ref_data, ref_count, out);

    simulation_result_free(&result);
}

static PassFlags evaluate_passes(const RunMetrics* base_metrics, real_t fw_start, const RunMetrics* dt_half_metrics,
    real_t rel_threshold, real_t night_day_ratio_max, real_t min_growth_rel,
    real_t dt_rgr_threshold_pct, real_t dt_biomass_threshold_pct,
    int use_external_data, int min_ref_points, real_t sucrose_nrmse_threshold, real_t starch_nrmse_threshold,
    real_t* out_dt_rgr_delta_pct, real_t* out_dt_biomass_delta_pct)
{
    PassFlags flags;
    memset(&flags, 0, sizeof(flags));

    flags.pass_rgr = base_metrics->rel_err <= rel_threshold;
    flags.pass_day_night = (base_metrics->day_avg_ph > REAL(0.0))
        && (base_metrics->night_avg_ph <= base_metrics->day_avg_ph * night_day_ratio_max);
    flags.pass_state = (base_metrics->nonfinite_samples == 0)
        && (base_metrics->min_sucrose >= REAL(0.0))
        && (base_metrics->min_starch >= REAL(0.0))
        && (base_metrics->min_total_biomass > REAL(0.0));
    flags.pass_growth = base_metrics->final_total_biomass >= fw_start * (REAL(1.0) + min_growth_rel);

    *out_dt_rgr_delta_pct = rel_pct(base_metrics->model_rgr, dt_half_metrics->model_rgr);
    *out_dt_biomass_delta_pct = rel_pct(base_metrics->final_total_biomass, dt_half_metrics->final_total_biomass);
    flags.pass_dt = (*out_dt_rgr_delta_pct <= dt_rgr_threshold_pct)
        && (*out_dt_biomass_delta_pct <= dt_biomass_threshold_pct);

    if (use_external_data) {
        flags.pass_external = (base_metrics->ref_points >= min_ref_points)
            && (base_metrics->sucrose_nrmse >= REAL(0.0))
            && (base_metrics->starch_nrmse >= REAL(0.0))
            && (base_metrics->sucrose_nrmse <= sucrose_nrmse_threshold)
            && (base_metrics->starch_nrmse <= starch_nrmse_threshold);
    } else {
        flags.pass_external = 1;
    }

    if (!flags.pass_rgr)
        append_fail_reason(flags.fail_reason, sizeof(flags.fail_reason), "rgr");
    if (!flags.pass_day_night)
        append_fail_reason(flags.fail_reason, sizeof(flags.fail_reason), "day_night_ph");
    if (!flags.pass_state)
        append_fail_reason(flags.fail_reason, sizeof(flags.fail_reason), "state_bounds");
    if (!flags.pass_growth)
        append_fail_reason(flags.fail_reason, sizeof(flags.fail_reason), "growth");
    if (!flags.pass_dt)
        append_fail_reason(flags.fail_reason, sizeof(flags.fail_reason), "dt_sensitivity");
    if (!flags.pass_external)
        append_fail_reason(flags.fail_reason, sizeof(flags.fail_reason), "sugar_starch_fit");

    flags.pass = flags.pass_rgr && flags.pass_day_night && flags.pass_state && flags.pass_growth && flags.pass_dt && flags.pass_external;
    if (flags.pass && flags.fail_reason[0] == '\0')
        snprintf(flags.fail_reason, sizeof(flags.fail_reason), "none");

    return flags;
}

int main(int argc, char** argv)
{
    const char* config_path = "config/default.config";
    const char* out_path = "validation_summary.csv";
    const char* sugars_starch_data_path = "data/paper/sugars_starch_reference.csv";

    int days = 30;
    real_t rel_threshold = REAL(0.10);
    real_t night_day_ratio_max = REAL(0.05);
    real_t min_growth_rel = REAL(0.02);
    real_t max_negative_step_abs = REAL(1e-10);
    real_t dt_rgr_threshold_pct = REAL(2.0);
    real_t dt_biomass_threshold_pct = REAL(3.0);
    real_t sucrose_nrmse_threshold = REAL(0.50);
    real_t starch_nrmse_threshold = REAL(0.25);
    int min_ref_points = 6;

    char overall_path[512];

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--config") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Missing value for --config\n");
                return 1;
            }
            config_path = argv[++i];
        } else if (strcmp(argv[i], "--days") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Missing value for --days\n");
                return 1;
            }
            days = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--out") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Missing value for --out\n");
                return 1;
            }
            out_path = argv[++i];
        } else if (strcmp(argv[i], "--threshold-rel") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Missing value for --threshold-rel\n");
                return 1;
            }
            rel_threshold = (real_t)atof(argv[++i]);
        } else if (strcmp(argv[i], "--threshold-night-day-ratio") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Missing value for --threshold-night-day-ratio\n");
                return 1;
            }
            night_day_ratio_max = (real_t)atof(argv[++i]);
        } else if (strcmp(argv[i], "--min-growth-rel") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Missing value for --min-growth-rel\n");
                return 1;
            }
            min_growth_rel = (real_t)atof(argv[++i]);
        } else if (strcmp(argv[i], "--max-negative-step-abs") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Missing value for --max-negative-step-abs\n");
                return 1;
            }
            max_negative_step_abs = (real_t)atof(argv[++i]);
        } else if (strcmp(argv[i], "--threshold-dt-rgr-pct") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Missing value for --threshold-dt-rgr-pct\n");
                return 1;
            }
            dt_rgr_threshold_pct = (real_t)atof(argv[++i]);
        } else if (strcmp(argv[i], "--threshold-dt-biomass-pct") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Missing value for --threshold-dt-biomass-pct\n");
                return 1;
            }
            dt_biomass_threshold_pct = (real_t)atof(argv[++i]);
        } else if (strcmp(argv[i], "--sugars-starch-data") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Missing value for --sugars-starch-data\n");
                return 1;
            }
            sugars_starch_data_path = argv[++i];
        } else if (strcmp(argv[i], "--threshold-sucrose-nrmse") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Missing value for --threshold-sucrose-nrmse\n");
                return 1;
            }
            sucrose_nrmse_threshold = (real_t)atof(argv[++i]);
        } else if (strcmp(argv[i], "--threshold-starch-nrmse") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Missing value for --threshold-starch-nrmse\n");
                return 1;
            }
            starch_nrmse_threshold = (real_t)atof(argv[++i]);
        } else if (strcmp(argv[i], "--min-ref-points") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Missing value for --min-ref-points\n");
                return 1;
            }
            min_ref_points = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
            usage(argv[0]);
            return 0;
        } else {
            fprintf(stderr, "Unknown argument: %s\n", argv[i]);
            usage(argv[0]);
            return 1;
        }
    }

    if (days < 2) {
        fprintf(stderr, "--days must be >= 2\n");
        return 1;
    }

    ReferencePoint ref_data[REF_POINTS_CAP];
    int ref_count = 0;
    int use_external_data = 0;

    if (load_reference_data(sugars_starch_data_path, ref_data, REF_POINTS_CAP, &ref_count) == 0 && ref_count > 0) {
        use_external_data = 1;
        printf("Loaded sugar/starch reference data: %d rows from %s\n", ref_count, sugars_starch_data_path);
    } else {
        printf("Warning: sugar/starch reference data unavailable at %s; external-fit checks disabled.\n", sugars_starch_data_path);
    }

    SimulationConfig base_cfg = simulation_default_config();
    Input base_input = generate_input();
    if (config_load_file(config_path, &base_cfg, &base_input) != 0) {
        fprintf(stderr, "Failed to load config: %s\n", config_path);
        return 1;
    }
    base_cfg.days = days;
    base_cfg.gui_enabled = 0;
    base_cfg.write_csv = 0;

    FILE* out = fopen(out_path, "w");
    if (!out) {
        fprintf(stderr, "Failed to open output CSV: %s\n", out_path);
        return 1;
    }
    fprintf(out,
        "photoperiod_h,days,model_rgr,lit_rgr,abs_error,rel_error_pct,pass_rgr,pass_day_night,pass_state,pass_growth,pass_dt,pass_external,pass,max_ph,day_avg_ph,night_avg_ph,min_sucrose,min_starch,min_total_biomass,max_biomass_drop,negative_biomass_steps,nonfinite_samples,final_total_biomass,dt_half_rgr,rgr_dt_delta_pct,dt_half_final_total_biomass,final_biomass_dt_delta_pct,ref_points,sucrose_scale,sucrose_nrmse,starch_scale,starch_nrmse,fail_reason\n");

    real_t sum_abs = REAL(0.0);
    real_t sum_sq = REAL(0.0);
    real_t sum_rel = REAL(0.0);
    real_t max_abs = REAL(0.0);

    int pass_count = 0;
    int pass_rgr_count = 0;
    int pass_day_night_count = 0;
    int pass_state_count = 0;
    int pass_growth_count = 0;
    int pass_dt_count = 0;
    int pass_external_count = 0;

    int n = (int)(sizeof(SCENARIOS) / sizeof(SCENARIOS[0]));

    for (int i = 0; i < n; ++i) {
        Input in = base_input;
        in.photo.photoperiod = (real_t)SCENARIOS[i].photoperiod_h;
        in.carbohydrates.lambda_sni = SCENARIOS[i].lambda_sni;
        in.growth.lambda_sb = SCENARIOS[i].lambda_sb;
        in.photo.feedback_on_photosynthesis = SCENARIOS[i].feedback;

        real_t fw_start = in.growth.total_biomass;

        RunMetrics base_metrics;
        memset(&base_metrics, 0, sizeof(base_metrics));
        compute_run_metrics(&base_cfg, &in, fw_start, SCENARIOS[i].lit_rgr,
            max_negative_step_abs, ref_data, ref_count, &base_metrics);

        SimulationConfig half_cfg = base_cfg;
        half_cfg.dt_hours = base_cfg.dt_hours * REAL(0.5);
        if (half_cfg.dt_hours <= REAL(1e-6))
            half_cfg.dt_hours = base_cfg.dt_hours;

        Input in_half = base_input;
        in_half.photo.photoperiod = (real_t)SCENARIOS[i].photoperiod_h;
        in_half.carbohydrates.lambda_sni = SCENARIOS[i].lambda_sni;
        in_half.growth.lambda_sb = SCENARIOS[i].lambda_sb;
        in_half.photo.feedback_on_photosynthesis = SCENARIOS[i].feedback;

        RunMetrics half_metrics;
        memset(&half_metrics, 0, sizeof(half_metrics));
        compute_run_metrics(&half_cfg, &in_half, fw_start, SCENARIOS[i].lit_rgr,
            max_negative_step_abs, NULL, 0, &half_metrics);

        real_t dt_rgr_delta_pct = REAL(0.0);
        real_t dt_biomass_delta_pct = REAL(0.0);
        PassFlags flags = evaluate_passes(&base_metrics, fw_start, &half_metrics,
            rel_threshold, night_day_ratio_max, min_growth_rel,
            dt_rgr_threshold_pct, dt_biomass_threshold_pct,
            use_external_data, min_ref_points, sucrose_nrmse_threshold, starch_nrmse_threshold,
            &dt_rgr_delta_pct, &dt_biomass_delta_pct);

        fprintf(out,
            "%d,%d,%.6f,%.6f,%.6f,%.3f,%d,%d,%d,%d,%d,%d,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6e,%d,%d,%.6f,%.6f,%.3f,%.6f,%.3f,%d,%.6f,%.6f,%.6f,%.6f,%s\n",
            SCENARIOS[i].photoperiod_h,
            base_cfg.days,
            (double)base_metrics.model_rgr,
            (double)SCENARIOS[i].lit_rgr,
            (double)base_metrics.abs_err,
            (double)(base_metrics.rel_err * REAL(100.0)),
            flags.pass_rgr,
            flags.pass_day_night,
            flags.pass_state,
            flags.pass_growth,
            flags.pass_dt,
            flags.pass_external,
            flags.pass,
            (double)base_metrics.max_ph,
            (double)base_metrics.day_avg_ph,
            (double)base_metrics.night_avg_ph,
            (double)base_metrics.min_sucrose,
            (double)base_metrics.min_starch,
            (double)base_metrics.min_total_biomass,
            (double)base_metrics.max_biomass_drop,
            base_metrics.negative_biomass_steps,
            base_metrics.nonfinite_samples,
            (double)base_metrics.final_total_biomass,
            (double)half_metrics.model_rgr,
            (double)dt_rgr_delta_pct,
            (double)half_metrics.final_total_biomass,
            (double)dt_biomass_delta_pct,
            base_metrics.ref_points,
            (double)base_metrics.sucrose_scale,
            (double)base_metrics.sucrose_nrmse,
            (double)base_metrics.starch_scale,
            (double)base_metrics.starch_nrmse,
            flags.fail_reason);

        sum_abs += base_metrics.abs_err;
        sum_sq += base_metrics.abs_err * base_metrics.abs_err;
        sum_rel += base_metrics.rel_err;
        if (base_metrics.abs_err > max_abs)
            max_abs = base_metrics.abs_err;

        if (flags.pass)
            pass_count++;
        if (flags.pass_rgr)
            pass_rgr_count++;
        if (flags.pass_day_night)
            pass_day_night_count++;
        if (flags.pass_state)
            pass_state_count++;
        if (flags.pass_growth)
            pass_growth_count++;
        if (flags.pass_dt)
            pass_dt_count++;
        if (flags.pass_external)
            pass_external_count++;
    }

    real_t mae = sum_abs / (real_t)n;
    real_t rmse = RSQRT(sum_sq / (real_t)n);
    real_t mape_pct = (sum_rel / (real_t)n) * REAL(100.0);
    real_t pass_rate_pct = ((real_t)pass_count / (real_t)n) * REAL(100.0);

    fprintf(out, "OVERALL,%d,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n", base_cfg.days);
    fclose(out);

    snprintf(overall_path, sizeof(overall_path), "%s", out_path);
    {
        char* ext = strrchr(overall_path, '.');
        if (ext && strcmp(ext, ".csv") == 0)
            *ext = '\0';
    }
    strncat(overall_path, "_overall.csv", sizeof(overall_path) - strlen(overall_path) - 1);

    FILE* out_overall = fopen(overall_path, "w");
    if (out_overall) {
        fprintf(out_overall,
            "days,mae,rmse,max_abs_error,mape_pct,pass_count,total_scenarios,pass_rate_pct,pass_rgr_count,pass_day_night_count,pass_state_count,pass_growth_count,pass_dt_count,pass_external_count,use_external_data,reference_rows,threshold_rel,threshold_pct,threshold_night_day_ratio_max,min_growth_rel_pct,threshold_dt_rgr_pct,threshold_dt_biomass_pct,max_negative_step_abs,threshold_sucrose_nrmse,threshold_starch_nrmse,min_ref_points\n");
        fprintf(out_overall,
            "%d,%.6f,%.6f,%.6f,%.3f,%d,%d,%.3f,%d,%d,%d,%d,%d,%d,%d,%d,%.6f,%.3f,%.6f,%.3f,%.3f,%.3f,%.6e,%.3f,%.3f,%d\n",
            base_cfg.days,
            (double)mae,
            (double)rmse,
            (double)max_abs,
            (double)mape_pct,
            pass_count,
            n,
            (double)pass_rate_pct,
            pass_rgr_count,
            pass_day_night_count,
            pass_state_count,
            pass_growth_count,
            pass_dt_count,
            pass_external_count,
            use_external_data,
            ref_count,
            (double)rel_threshold,
            (double)(rel_threshold * REAL(100.0)),
            (double)night_day_ratio_max,
            (double)(min_growth_rel * REAL(100.0)),
            (double)dt_rgr_threshold_pct,
            (double)dt_biomass_threshold_pct,
            (double)max_negative_step_abs,
            (double)sucrose_nrmse_threshold,
            (double)starch_nrmse_threshold,
            min_ref_points);
        fclose(out_overall);
    }

    printf("Validation written to %s\n", out_path);
    printf("Overall metrics written to %s\n", overall_path);
    printf("MAE=%.6f RMSE=%.6f MAPE=%.3f%% pass=%d/%d (RGR threshold=%.3f%%, dt RGR threshold=%.3f%%, dt biomass threshold=%.3f%%, sucrose nRMSE threshold=%.3f, starch nRMSE threshold=%.3f)\n",
        (double)mae,
        (double)rmse,
        (double)mape_pct,
        pass_count,
        n,
        (double)(rel_threshold * REAL(100.0)),
        (double)dt_rgr_threshold_pct,
        (double)dt_biomass_threshold_pct,
        (double)sucrose_nrmse_threshold,
        (double)starch_nrmse_threshold);
    return 0;
}

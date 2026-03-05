#define SOLVER_REAL_AS double
#include "config/config.h"
#include "model/input.h"
#include "simulation/simulation.h"

#include <errno.h>
#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>

typedef struct {
    int photoperiod_h;
    real_t lit_rgr;
} Scenario;

typedef struct {
    char name[96];
    int enabled;
    real_t min_value;
    real_t max_value;
    int log_scale;
    size_t offset;
    real_t base_value;
} SweepParam;

typedef struct {
    int is_set;
    char param[96];
    real_t value;
    int photoperiod_h;
    real_t model_rgr;
    real_t lit_rgr;
    real_t abs_err;
} BestByPhotoperiod;

typedef struct {
    int ok;
    real_t rgr;
} RunResultMsg;

static const Scenario DEFAULT_SCENARIOS[] = {
    { 4, REAL(0.0680) },
    { 6, REAL(0.1135) },
    { 8, REAL(0.1708) },
    { 12, REAL(0.2600) },
    // { 18, REAL(0.3065) }
};

static char* trim(char* s)
{
    while (*s == ' ' || *s == '\t' || *s == '\n' || *s == '\r')
        s++;
    if (*s == '\0')
        return s;
    char* end = s + strlen(s) - 1;
    while (end > s && (*end == ' ' || *end == '\t' || *end == '\n' || *end == '\r')) {
        *end = '\0';
        end--;
    }
    return s;
}

static int parse_real(const char* s, real_t* out)
{
    if (!s || !out)
        return -1;
    errno = 0;
    char* end = NULL;
    double v = strtod(s, &end);
    if (errno != 0 || end == s)
        return -1;
    while (*end == ' ' || *end == '\t')
        end++;
    if (*end != '\0')
        return -1;
    *out = (real_t)v;
    return 0;
}

static int parse_int(const char* s, int* out)
{
    if (!s || !out)
        return -1;
    errno = 0;
    char* end = NULL;
    long v = strtol(s, &end, 10);
    if (errno != 0 || end == s)
        return -1;
    while (*end == ' ' || *end == '\t')
        end++;
    if (*end != '\0')
        return -1;
    *out = (int)v;
    return 0;
}

static int is_log_scale(const char* s)
{
    return strcmp(s, "log") == 0 || strcmp(s, "LOG") == 0 || strcmp(s, "Log") == 0;
}

static int load_params_csv(
    const char* path,
    const Input* base_input,
    SweepParam* params,
    int max_params)
{
    FILE* f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "Could not open params CSV: %s\n", path);
        return -1;
    }

    int count = 0;
    char line[512];
    int line_no = 0;
    while (fgets(line, sizeof(line), f)) {
        line_no++;
        char* s = trim(line);
        if (*s == '\0' || *s == '#')
            continue;

        if (line_no == 1 && strstr(s, "name,enabled,min,max,scale") == s)
            continue;

        char* cols[5] = { 0 };
        int col_count = 0;
        char* tok = strtok(s, ",");
        while (tok && col_count < 5) {
            cols[col_count++] = trim(tok);
            tok = strtok(NULL, ",");
        }
        if (col_count != 5) {
            fprintf(stderr, "Invalid params CSV row at line %d\n", line_no);
            fclose(f);
            return -1;
        }

        if (count >= max_params) {
            fprintf(stderr, "Too many params in CSV (max=%d)\n", max_params);
            fclose(f);
            return -1;
        }

        SweepParam p;
        memset(&p, 0, sizeof(p));
        snprintf(p.name, sizeof(p.name), "%s", cols[0]);
        if (parse_int(cols[1], &p.enabled) != 0 || parse_real(cols[2], &p.min_value) != 0 || parse_real(cols[3], &p.max_value) != 0) {
            fprintf(stderr, "Invalid numeric values at line %d\n", line_no);
            fclose(f);
            return -1;
        }
        p.log_scale = is_log_scale(cols[4]);

        if (config_get_input_field_offset(p.name, &p.offset) != 0) {
            fprintf(stderr, "Unknown input field in sweep CSV: '%s' (line %d)\n", p.name, line_no);
            fclose(f);
            return -1;
        }

        p.base_value = *(const real_t*)((const char*)base_input + p.offset);
        params[count++] = p;
    }

    fclose(f);
    return count;
}

static real_t grid_value(const SweepParam* p, int index, int points)
{
    if (points <= 1)
        return p->min_value;

    real_t alpha = (real_t)index / (real_t)(points - 1);
    if (!p->log_scale) {
        return p->min_value + alpha * (p->max_value - p->min_value);
    }

    real_t lo = RLOG(RMAX(p->min_value, REAL(1e-12)));
    real_t hi = RLOG(RMAX(p->max_value, REAL(1e-12)));
    return REXP(lo + alpha * (hi - lo));
}

static void print_usage(const char* argv0)
{
    printf("Usage: %s [--config path] [--params csv] [--out csv] [--days N] [--points N] [--dt-hours H] [--solver-tol X] [--run-timeout-sec N]\n", argv0);
}

static void print_progress(int done, int total, time_t start_ts)
{
    if (total <= 0)
        return;

    const int bar_w = 36;
    double frac = (double)done / (double)total;
    if (frac < 0.0)
        frac = 0.0;
    if (frac > 1.0)
        frac = 1.0;

    int filled = (int)(frac * bar_w);
    time_t now = time(NULL);
    double elapsed = difftime(now, start_ts);
    double eta = 0.0;
    if (done > 0)
        eta = (elapsed / (double)done) * (double)(total - done);

    int eta_h = (int)(eta / 3600.0);
    int eta_m = (int)((eta - eta_h * 3600) / 60.0);
    int eta_s = (int)(eta - eta_h * 3600 - eta_m * 60);

    printf("\r[");
    for (int i = 0; i < bar_w; ++i)
        putchar(i < filled ? '#' : '-');
    printf("] %6.2f%%  %d/%d  ETA %02d:%02d:%02d",
        frac * 100.0, done, total, eta_h, eta_m, eta_s);
    fflush(stdout);
}

static void build_best_output_path(const char* out_path, char* dst, size_t dst_size)
{
    const char* ext = strrchr(out_path, '.');
    if (ext && strcmp(ext, ".csv") == 0) {
        size_t prefix_len = (size_t)(ext - out_path);
        if (prefix_len + strlen("_best_by_photoperiod.csv") + 1 <= dst_size) {
            memcpy(dst, out_path, prefix_len);
            dst[prefix_len] = '\0';
            strcat(dst, "_best_by_photoperiod.csv");
            return;
        }
    }
    snprintf(dst, dst_size, "%s_best_by_photoperiod.csv", out_path);
}

static int run_simulation_with_timeout(
    const SimulationConfig* cfg,
    const Input* input,
    int timeout_sec,
    real_t* rgr_out)
{
    if (!cfg || !input || !rgr_out)
        return -1;

    int fds[2];
    if (pipe(fds) != 0)
        return -1;

    pid_t pid = fork();
    if (pid < 0) {
        close(fds[0]);
        close(fds[1]);
        return -1;
    }

    if (pid == 0) {
        close(fds[0]);

        RunResultMsg msg;
        msg.ok = 0;
        msg.rgr = REAL(0.0);

        Input run_input = *input;
        real_t fw_start = run_input.core.total_biomass;
        SimulationResult result;
        simulate_days(cfg, &run_input, &result);
        msg.rgr = compute_rgr(fw_start, run_input.core.total_biomass, (real_t)cfg->days);
        msg.ok = 1;
        simulation_result_free(&result);

        (void)write(fds[1], &msg, sizeof(msg));
        close(fds[1]);
        _exit(0);
    }

    close(fds[1]);

    time_t start = time(NULL);
    int status = 0;
    int finished = 0;
    while (!finished) {
        pid_t w = waitpid(pid, &status, WNOHANG);
        if (w == pid) {
            finished = 1;
            break;
        }
        if (w < 0) {
            kill(pid, SIGKILL);
            waitpid(pid, NULL, 0);
            close(fds[0]);
            return -1;
        }

        if (difftime(time(NULL), start) >= timeout_sec) {
            kill(pid, SIGKILL);
            waitpid(pid, NULL, 0);
            close(fds[0]);
            return 1; // timeout
        }

        usleep(100000);
    }

    RunResultMsg msg;
    ssize_t nread = read(fds[0], &msg, sizeof(msg));
    close(fds[0]);

    if (!WIFEXITED(status) || WEXITSTATUS(status) != 0)
        return -1;
    if (nread != (ssize_t)sizeof(msg) || !msg.ok)
        return -1;

    *rgr_out = msg.rgr;
    return 0;
}

int main(int argc, char** argv)
{
    const char* config_path = "config/default.config";
    const char* params_path = "config/sweep_params.csv";
    const char* out_path = "sweep_results.csv";
    int points = 9;
    int days_override = 30;
    real_t dt_hours_override = REAL(0.05);
    real_t solver_tol_override = REAL(1e-6);
    int run_timeout_sec = 30;

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--config") == 0 && i + 1 < argc) {
            config_path = argv[++i];
        } else if (strcmp(argv[i], "--params") == 0 && i + 1 < argc) {
            params_path = argv[++i];
        } else if (strcmp(argv[i], "--out") == 0 && i + 1 < argc) {
            out_path = argv[++i];
        } else if (strcmp(argv[i], "--points") == 0 && i + 1 < argc) {
            points = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--days") == 0 && i + 1 < argc) {
            days_override = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--dt-hours") == 0 && i + 1 < argc) {
            dt_hours_override = (real_t)atof(argv[++i]);
        } else if (strcmp(argv[i], "--solver-tol") == 0 && i + 1 < argc) {
            solver_tol_override = (real_t)atof(argv[++i]);
        } else if (strcmp(argv[i], "--run-timeout-sec") == 0 && i + 1 < argc) {
            run_timeout_sec = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
            print_usage(argv[0]);
            return 0;
        } else {
            fprintf(stderr, "Unknown argument: %s\n", argv[i]);
            print_usage(argv[0]);
            return 1;
        }
    }

    if (points < 2)
        points = 2;
    if (days_override <= 0)
        days_override = 1;
    if (dt_hours_override <= REAL(0.0))
        dt_hours_override = REAL(0.05);
    if (solver_tol_override <= REAL(0.0))
        solver_tol_override = REAL(1e-6);
    if (run_timeout_sec <= 0)
        run_timeout_sec = 30;

    SimulationConfig base_cfg = simulation_default_config();
    Input base_input = generate_input();
    if (config_load_file(config_path, &base_cfg, &base_input) != 0)
        return 1;

    base_cfg.days = days_override;
    base_cfg.dt_hours = dt_hours_override;
    base_cfg.solver_tolerance = solver_tol_override;
    base_cfg.gui_enabled = 0;
    base_cfg.write_csv = 0;

    SweepParam params[128];
    int param_count = load_params_csv(params_path, &base_input, params, 128);
    if (param_count <= 0) {
        fprintf(stderr, "No sweep params loaded from: %s\n", params_path);
        return 1;
    }

    FILE* out = fopen(out_path, "w");
    if (!out) {
        fprintf(stderr, "Could not open output CSV: %s\n", out_path);
        return 1;
    }

    fprintf(out, "param,value,scenario_photoperiod,model_rgr,lit_rgr,abs_error,days\n");

    const int scenario_count = (int)(sizeof(DEFAULT_SCENARIOS) / sizeof(DEFAULT_SCENARIOS[0]));
    int enabled_count = 0;
    for (int i = 0; i < param_count; ++i) {
        if (params[i].enabled)
            enabled_count++;
    }
    int total_work = enabled_count * points * scenario_count;
    int done_work = 0;
    int total_runs = 0;
    time_t start_ts = time(NULL);

    BestByPhotoperiod best_per_scenario[sizeof(DEFAULT_SCENARIOS) / sizeof(DEFAULT_SCENARIOS[0])] = { 0 };
    printf("Sweep config: enabled_params=%d points=%d scenarios=%d total_runs=%d days=%d dt_hours=%g solver_tol=%g timeout=%ds\n",
        enabled_count, points, scenario_count, total_work, base_cfg.days, (double)base_cfg.dt_hours, (double)base_cfg.solver_tolerance, run_timeout_sec);
    print_progress(0, total_work, start_ts);

    for (int p = 0; p < param_count; ++p) {
        if (!params[p].enabled)
            continue;

        real_t best_mae = REAL(1e9);
        real_t best_value = params[p].base_value;

        for (int vi = 0; vi < points; ++vi) {
            real_t val = grid_value(&params[p], vi, points);
            real_t sum_abs_err = REAL(0.0);

            for (int si = 0; si < scenario_count; ++si) {
                Input run_input = base_input;
                SimulationConfig run_cfg = base_cfg;
                run_input.core.photoperiod = (real_t)DEFAULT_SCENARIOS[si].photoperiod_h;

                real_t* dst = (real_t*)((char*)&run_input + params[p].offset);
                *dst = val;

                real_t rgr = REAL(0.0);
                int run_rc = run_simulation_with_timeout(&run_cfg, &run_input, run_timeout_sec, &rgr);
                int timed_out = (run_rc == 1);
                int failed = (run_rc != 0);
                real_t abs_err = failed ? REAL(1e6) : RABS(rgr - DEFAULT_SCENARIOS[si].lit_rgr);
                sum_abs_err += abs_err;

                if (!best_per_scenario[si].is_set || abs_err < best_per_scenario[si].abs_err) {
                    best_per_scenario[si].is_set = 1;
                    size_t param_len = strlen(params[p].name);
                    if (param_len >= sizeof(best_per_scenario[si].param))
                        param_len = sizeof(best_per_scenario[si].param) - 1;
                    memcpy(best_per_scenario[si].param, params[p].name, param_len);
                    best_per_scenario[si].param[param_len] = '\0';
                    best_per_scenario[si].value = val;
                    best_per_scenario[si].photoperiod_h = DEFAULT_SCENARIOS[si].photoperiod_h;
                    best_per_scenario[si].model_rgr = rgr;
                    best_per_scenario[si].lit_rgr = DEFAULT_SCENARIOS[si].lit_rgr;
                    best_per_scenario[si].abs_err = abs_err;
                }

                fprintf(out, "%s,%.10g,%d,%.10g,%.10g,%.10g,%d\n",
                    params[p].name,
                    (double)val,
                    DEFAULT_SCENARIOS[si].photoperiod_h,
                    (double)rgr,
                    (double)DEFAULT_SCENARIOS[si].lit_rgr,
                    (double)abs_err,
                    run_cfg.days);
                if (timed_out) {
                    fprintf(stderr, "\nWarning: timeout for param=%s value=%g photoperiod=%dh (>%ds)\n",
                        params[p].name, (double)val, DEFAULT_SCENARIOS[si].photoperiod_h, run_timeout_sec);
                } else if (failed) {
                    fprintf(stderr, "\nWarning: failed run for param=%s value=%g photoperiod=%dh\n",
                        params[p].name, (double)val, DEFAULT_SCENARIOS[si].photoperiod_h);
                }
                total_runs++;
                done_work++;
                print_progress(done_work, total_work, start_ts);
            }

            real_t mae = sum_abs_err / (real_t)(sizeof(DEFAULT_SCENARIOS) / sizeof(DEFAULT_SCENARIOS[0]));
            if (mae < best_mae) {
                best_mae = mae;
                best_value = val;
            }
        }

        printf("best %-45s value=%g  mae=%g\n", params[p].name, (double)best_value, (double)best_mae);
    }
    printf("\n");

    fclose(out);

    char best_out_path[512];
    build_best_output_path(out_path, best_out_path, sizeof(best_out_path));
    FILE* best_out = fopen(best_out_path, "w");
    if (!best_out) {
        fprintf(stderr, "Could not open best-per-photoperiod CSV: %s\n", best_out_path);
        return 1;
    }
    fprintf(best_out, "scenario_photoperiod,best_param,best_value,model_rgr,lit_rgr,abs_error\n");
    printf("Best by photoperiod:\n");
    for (int si = 0; si < scenario_count; ++si) {
        if (!best_per_scenario[si].is_set)
            continue;
        fprintf(best_out, "%d,%s,%.10g,%.10g,%.10g,%.10g\n",
            best_per_scenario[si].photoperiod_h,
            best_per_scenario[si].param,
            (double)best_per_scenario[si].value,
            (double)best_per_scenario[si].model_rgr,
            (double)best_per_scenario[si].lit_rgr,
            (double)best_per_scenario[si].abs_err);
        printf("  %2dh: %-45s value=%g abs_err=%g\n",
            best_per_scenario[si].photoperiod_h,
            best_per_scenario[si].param,
            (double)best_per_scenario[si].value,
            (double)best_per_scenario[si].abs_err);
    }
    fclose(best_out);

    printf("Sweep finished. runs=%d output=%s best=%s\n", total_runs, out_path, best_out_path);
    return 0;
}

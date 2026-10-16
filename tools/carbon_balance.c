/*
 * carbon_balance.c — closure diagnostics for the carbohydrate subsystem
 *
 * Purpose
 * -------
 * Reviewers asked whether the carbon budget of the model closes (issues #16,
 * #17).  The question cannot be answered by post-processing the ODE output:
 * the state trajectory records pool *sizes*, whereas closure is a statement
 * about time-integrated *fluxes*.  This tool therefore integrates the
 * carbohydrate fluxes alongside the state, using the same RK4 weights and the
 * same step size, so that the resulting budget is exact to the arithmetic of
 * the integrator rather than to the accuracy of a quadrature applied
 * afterwards.
 *
 * Scope of the claim
 * ------------------
 * Closure is asserted for the carbohydrate subsystem only:
 *
 *     {a (starch), s (cytosolic sucrose), s_v (vacuolar sucrose)}
 *
 * The biomass equations are deliberately excluded.  In the model
 *
 *     dB_l/dt = lambda_sb (1 - phi_r) phi_load E_night B_l - ...
 *
 * lambda_sb carries units of gFW^-1, so the growth terms are a *rate law*
 * calibrated against observed RGR, not a stoichiometric conversion of sucrose
 * into biomass carbon.  The sucrose drain phi_load * E_night and the biomass
 * gain lambda_sb * phi_load * E_night * B_l differ by the factor
 * lambda_sb * B_l, which equals unity only at B_l = 1/lambda_sb ~ 170 gFW,
 * four orders of magnitude above the seedling biomass simulated here.  A
 * whole-plant carbon budget would require a stoichiometric growth equation;
 * see docs/revision/math_audit.md section 5.  What this tool verifies is that
 * within its stated scope no carbon is created or destroyed numerically.
 *
 * Budget identity (per pool, over [0, T])
 * ---------------------------------------
 *   starch:   a(T) - a(0)   = INT k_s * A_n            - INT D
 *   sucrose:  s(T) - s(0)   = INT (1-k_s) * A_n + INT D
 *                             - INT C_up - INT C_tr - INT r_m s
 *                             - INT phi_load E_night
 *                             - INT k_in s + INT k_out s_v
 *   vacuole:  s_v(T)-s_v(0) = INT k_in s - INT k_out s_v
 *
 * and for the subsystem total (the vacuole exchange cancels identically):
 *   (a+s+s_v)(T) - (a+s+s_v)(0) = INT A_n - INT (C_up + C_tr + r_m s + phi_load E_night)
 *
 * Residual = LHS - RHS.  Reported both absolute and relative to gross
 * assimilation INT A_n.
 *
 * Output
 * ------
 *   results/carbon_balance/carbon_balance.csv       per-photoperiod totals
 *   results/carbon_balance/carbon_daily_<PP>h.csv   per-day fluxes
 */

#define _POSIX_C_SOURCE 200809L

#include "config/config.h"
#include "model/input.h"
#include "model/model.h"
#include "model/state.h"
#include "simulation/dynamics.h"
#include "simulation/simulation.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

/* ------------------------------------------------------------------ flux accumulator layout */

typedef enum {
    F_ASSIM = 0,      /* A_n                      gross assimilation        */
    F_TO_STARCH,      /* k_s A_n                  assimilate -> starch      */
    F_TO_SUCROSE,     /* (1-k_s) A_n              assimilate -> sucrose     */
    F_DEGRADE,        /* D                        starch -> sucrose         */
    F_UPTAKE,         /* C_up                     nutrient uptake cost      */
    F_TRANSPORT,      /* C_tr                     phloem transport cost     */
    F_MAINT,          /* r_m s                    maintenance respiration   */
    F_LOAD,           /* phi_load E_night         export to growth sink     */
    F_VAC_IN,         /* k_in s                   cytosol -> vacuole        */
    F_VAC_OUT,        /* k_out s_v                vacuole -> cytosol        */
    F_DIM
} FluxIndex;

static const char* FLUX_NAME[F_DIM] = {
    "assimilation", "to_starch", "to_sucrose", "starch_degradation",
    "uptake_cost", "transport_cost", "maintenance_resp", "phloem_load",
    "vacuole_in", "vacuole_out"
};

#define AUG_DIM (X_DIM + F_DIM)

/* ------------------------------------------------------------------ augmented RHS */

/*
 * Evaluates the model RHS on the first X_DIM components and, from the same
 * algebraic state, the instantaneous carbohydrate fluxes on the remaining
 * F_DIM components.  Because both are produced by one call, the flux
 * accumulators see exactly the same stage values as the state integration.
 */
static void augmented_rhs(real_t t, const real_t* x, real_t* dx, PlantCtx* ctx)
{
    simulation_plant_rhs(t, x, dx, ctx);

    /* Recreate the algebraic state at this stage to read the fluxes. */
    Input s = *(ctx->base);
    simulation_unpack_state(&s, x);
    s.growth.total_biomass = s.growth.leaf_biomass + s.growth.root_biomass;
    simulation_compute_algebraic(&s, fmod(t, REAL(24.0)), ctx->starch_night_start);

    real_t A_n   = s.photo.photosynthesis;
    real_t k_s   = s.carbohydrates.starch_partition_coeff;
    real_t D     = s.carbohydrates.starch_degradation_rate;
    real_t suc   = s.carbohydrates.sucrose;
    real_t vac   = s.carbohydrates.sucrose_vacuole;

    real_t* f = dx + X_DIM;
    f[F_ASSIM]      = A_n;
    f[F_TO_STARCH]  = k_s * A_n;
    f[F_TO_SUCROSE] = (REAL(1.0) - k_s) * A_n;
    f[F_DEGRADE]    = D;
    f[F_UPTAKE]     = s.carbohydrates.uptake_cost;
    f[F_TRANSPORT]  = s.carbohydrates.transport_cost;
    f[F_MAINT]      = s.carbohydrates.respiration_frequency * suc;
    f[F_LOAD]       = s.carbohydrates.sucrose_loading_frequency
                    * s.carbohydrates.night_efficiency_starch;
    f[F_VAC_IN]     = s.carbohydrates.k_vacuole_in  * suc;
    f[F_VAC_OUT]    = s.carbohydrates.k_vacuole_out * vac;
}

/* Classical RK4 step on the augmented system. */
static void rk4_step(real_t t, real_t h, real_t* x, PlantCtx* ctx,
                     real_t* k1, real_t* k2, real_t* k3, real_t* k4, real_t* tmp)
{
    augmented_rhs(t, x, k1, ctx);

    for (size_t i = 0; i < AUG_DIM; ++i) tmp[i] = x[i] + REAL(0.5) * h * k1[i];
    augmented_rhs(t + REAL(0.5) * h, tmp, k2, ctx);

    for (size_t i = 0; i < AUG_DIM; ++i) tmp[i] = x[i] + REAL(0.5) * h * k2[i];
    augmented_rhs(t + REAL(0.5) * h, tmp, k3, ctx);

    for (size_t i = 0; i < AUG_DIM; ++i) tmp[i] = x[i] + h * k3[i];
    augmented_rhs(t + h, tmp, k4, ctx);

    for (size_t i = 0; i < AUG_DIM; ++i)
        x[i] += (h / REAL(6.0)) * (k1[i] + REAL(2.0) * k2[i] + REAL(2.0) * k3[i] + k4[i]);
}

/* ------------------------------------------------------------------ per-photoperiod run */

typedef struct {
    int    photoperiod_h;
    double flux[F_DIM];
    double a0, a1, s0, s1, v0, v1;
    double res_starch, res_sucrose, res_vacuole, res_total;
    double rel_total;
    double bl0, bl1, br0, br1;
} BalanceRow;

static void ensure_dir(const char* path)
{
    mkdir(path, 0775);
}

static int run_photoperiod(const SimulationConfig* cfg_in, const Input* base,
                           int ph_h, real_t dt, BalanceRow* row,
                           const char* daily_path)
{
    Input input = *base;
    input.photo.photoperiod = (real_t)ph_h;

    /* Same soil-dependent rescaling of lambda_sb that simulate_days applies,
     * so this tool integrates exactly the production model. */
    input.growth.lambda_sb = lambda_sb_f(
        input.growth.lambda_sb,
        input.nutrients.nitrogen_soil_content,
        input.nutrients.phosphorus_soil_content);

    PlantCtx ctx = {
        .base = &input,
        .starch_night_start = input.carbohydrates.starch,
        .was_light = 0
    };

    Input probe = input;
    update_light_conditions(&probe, REAL(0.0));
    ctx.was_light = (probe.photo.light != REAL(0.0));

    real_t x[AUG_DIM];
    memset(x, 0, sizeof(x));
    {
        PlantState st = plant_state_from_input(&input);
        plant_state_pack(x, &st);
    }

    row->photoperiod_h = ph_h;
    row->a0 = (double)x[X_STARCH];
    row->s0 = (double)x[X_SUCROSE];
    row->v0 = (double)x[X_SUCROSE_VAC];
    row->bl0 = (double)x[X_LEAF_BIOMASS];
    row->br0 = (double)x[X_ROOT_BIOMASS];

    real_t k1[AUG_DIM], k2[AUG_DIM], k3[AUG_DIM], k4[AUG_DIM], tmp[AUG_DIM];

    FILE* fd = daily_path ? fopen(daily_path, "w") : NULL;
    if (fd) {
        fprintf(fd, "day,starch,sucrose,sucrose_vacuole,leaf_biomass,root_biomass");
        for (int i = 0; i < F_DIM; ++i) fprintf(fd, ",d_%s", FLUX_NAME[i]);
        fprintf(fd, ",residual_total,residual_rel\n");
    }

    int    days       = cfg_in->days;
    long   steps_day  = (long)llround(24.0 / (double)dt);
    double prev_flux[F_DIM] = {0};
    double prev_pool  = (double)(x[X_STARCH] + x[X_SUCROSE] + x[X_SUCROSE_VAC]);

    for (int d = 0; d < days; ++d) {
        for (long k = 0; k < steps_day; ++k) {
            real_t t = (real_t)((double)d * 24.0 + (double)k * (double)dt);

            /* Night-start starch reference: the production observer latches the
             * starch level at the light -> dark transition. */
            Input probe2 = input;
            update_light_conditions(&probe2, fmod((double)t, 24.0));
            int is_light = (probe2.photo.light != REAL(0.0));
            if (ctx.was_light && !is_light)
                ctx.starch_night_start = x[X_STARCH];
            ctx.was_light = is_light;

            rk4_step(t, dt, x, &ctx, k1, k2, k3, k4, tmp);
        }

        if (fd) {
            double pool = (double)(x[X_STARCH] + x[X_SUCROSE] + x[X_SUCROSE_VAC]);
            double gain = (double)x[X_DIM + F_ASSIM] - prev_flux[F_ASSIM];
            double loss = ((double)x[X_DIM + F_UPTAKE]    - prev_flux[F_UPTAKE])
                        + ((double)x[X_DIM + F_TRANSPORT] - prev_flux[F_TRANSPORT])
                        + ((double)x[X_DIM + F_MAINT]     - prev_flux[F_MAINT])
                        + ((double)x[X_DIM + F_LOAD]      - prev_flux[F_LOAD]);
            double res = (pool - prev_pool) - (gain - loss);

            fprintf(fd, "%d,%.10g,%.10g,%.10g,%.10g,%.10g",
                d + 1,
                (double)x[X_STARCH], (double)x[X_SUCROSE], (double)x[X_SUCROSE_VAC],
                (double)x[X_LEAF_BIOMASS], (double)x[X_ROOT_BIOMASS]);
            for (int i = 0; i < F_DIM; ++i) {
                fprintf(fd, ",%.10g", (double)x[X_DIM + i] - prev_flux[i]);
                prev_flux[i] = (double)x[X_DIM + i];
            }
            fprintf(fd, ",%.6e,%.6e\n", res, (gain > 0.0) ? res / gain : 0.0);
            prev_pool = pool;
        }
    }

    if (fd) fclose(fd);

    for (int i = 0; i < F_DIM; ++i) row->flux[i] = (double)x[X_DIM + i];
    row->a1 = (double)x[X_STARCH];
    row->s1 = (double)x[X_SUCROSE];
    row->v1 = (double)x[X_SUCROSE_VAC];
    row->bl1 = (double)x[X_LEAF_BIOMASS];
    row->br1 = (double)x[X_ROOT_BIOMASS];

    /* Per-pool residuals. */
    row->res_starch = (row->a1 - row->a0)
                    - (row->flux[F_TO_STARCH] - row->flux[F_DEGRADE]);

    row->res_sucrose = (row->s1 - row->s0)
                     - (row->flux[F_TO_SUCROSE] + row->flux[F_DEGRADE]
                        - row->flux[F_UPTAKE] - row->flux[F_TRANSPORT]
                        - row->flux[F_MAINT]  - row->flux[F_LOAD]
                        - row->flux[F_VAC_IN] + row->flux[F_VAC_OUT]);

    row->res_vacuole = (row->v1 - row->v0)
                     - (row->flux[F_VAC_IN] - row->flux[F_VAC_OUT]);

    double pool_change = (row->a1 + row->s1 + row->v1)
                       - (row->a0 + row->s0 + row->v0);
    double net = row->flux[F_ASSIM]
               - (row->flux[F_UPTAKE] + row->flux[F_TRANSPORT]
                  + row->flux[F_MAINT] + row->flux[F_LOAD]);
    row->res_total = pool_change - net;
    row->rel_total = (row->flux[F_ASSIM] > 0.0)
                   ? row->res_total / row->flux[F_ASSIM] : 0.0;

    return 0;
}

/* ------------------------------------------------------------------ main */

int main(int argc, char** argv)
{
    const char* config_path = "config/default.config";
    const char* out_dir     = "results/carbon_balance";
    int    days_override    = 0;
    double dt_hours         = 0.0;   /* 0 => take from config */

    for (int i = 1; i < argc; ++i) {
        if      (strcmp(argv[i], "--config") == 0 && i + 1 < argc) config_path = argv[++i];
        else if (strcmp(argv[i], "--out")    == 0 && i + 1 < argc) out_dir     = argv[++i];
        else if (strcmp(argv[i], "--days")   == 0 && i + 1 < argc) days_override = atoi(argv[++i]);
        else if (strcmp(argv[i], "--dt")     == 0 && i + 1 < argc) dt_hours    = atof(argv[++i]);
        else {
            fprintf(stderr,
                "Usage: %s [--config <path>] [--out <dir>] [--days N] [--dt hours]\n", argv[0]);
            return 1;
        }
    }

    /* config_load_file() assigns only the keys it finds in the file, so both
     * structs must carry their defaults first.  Without this, every key absent
     * from config/default.config -- including carbohydrates.sucrose_vacuole and
     * k_vacuole_in/out -- retains uninitialised stack contents. */
    SimulationConfig cfg   = simulation_default_config();
    Input            input = generate_input();
    if (config_load_file(config_path, &cfg, &input) != 0) {
        fprintf(stderr, "Failed to load config: %s\n", config_path);
        return 1;
    }
    cfg.gui_enabled = 0;
    cfg.write_csv   = 0;
    if (days_override > 0) cfg.days = days_override;

    real_t dt = (dt_hours > 0.0) ? (real_t)dt_hours : cfg.dt_hours;
    if (!(dt > REAL(0.0))) dt = REAL(0.01);

    ensure_dir("results");
    ensure_dir(out_dir);

    printf("=== Carbohydrate carbon-balance closure ===\n");
    printf("  config: %s\n  days: %d\n  dt: %g h  (%ld steps/day)\n",
        config_path, cfg.days, (double)dt, (long)llround(24.0 / (double)dt));
    printf("  scope: {starch, sucrose_cytosol, sucrose_vacuole}\n");
    printf("  biomass excluded by construction: lambda_sb has units gFW^-1\n\n");

    const int PPS[] = { 4, 6, 8, 12 };
    const int NPP   = (int)(sizeof(PPS) / sizeof(PPS[0]));
    BalanceRow rows[4];

    for (int i = 0; i < NPP; ++i) {
        char daily[512];
        snprintf(daily, sizeof(daily), "%s/carbon_daily_%dh.csv", out_dir, PPS[i]);
        if (run_photoperiod(&cfg, &input, PPS[i], dt, &rows[i], daily) != 0) {
            fprintf(stderr, "run failed for %d h\n", PPS[i]);
            return 1;
        }
        printf("  %2dh  INT A_n = %12.5f   residual = %11.4e   relative = %11.4e\n",
            PPS[i], rows[i].flux[F_ASSIM], rows[i].res_total, rows[i].rel_total);
    }

    /* ---- summary CSV ---- */
    char path[512];
    snprintf(path, sizeof(path), "%s/carbon_balance.csv", out_dir);
    FILE* f = fopen(path, "w");
    if (!f) { fprintf(stderr, "cannot write %s\n", path); return 1; }

    fprintf(f, "photoperiod_h,days,dt_hours");
    for (int i = 0; i < F_DIM; ++i) fprintf(f, ",int_%s", FLUX_NAME[i]);
    fprintf(f, ",starch_start,starch_end,sucrose_start,sucrose_end,"
               "vacuole_start,vacuole_end,leaf_biomass_start,leaf_biomass_end,"
               "root_biomass_start,root_biomass_end,"
               "residual_starch,residual_sucrose,residual_vacuole,"
               "residual_total,residual_relative\n");

    for (int i = 0; i < NPP; ++i) {
        BalanceRow* r = &rows[i];
        fprintf(f, "%d,%d,%.10g", r->photoperiod_h, cfg.days, (double)dt);
        for (int k = 0; k < F_DIM; ++k) fprintf(f, ",%.12g", r->flux[k]);
        fprintf(f, ",%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g",
            r->a0, r->a1, r->s0, r->s1, r->v0, r->v1,
            r->bl0, r->bl1, r->br0, r->br1);
        fprintf(f, ",%.6e,%.6e,%.6e,%.6e,%.6e\n",
            r->res_starch, r->res_sucrose, r->res_vacuole,
            r->res_total, r->rel_total);
    }
    fclose(f);
    printf("\n  -> %s\n", path);

    /* ---- the scope caveat, quantified ----
     * Prints the factor by which the biomass gain term differs from the
     * sucrose drain it is nominally paired with, so the manuscript can state
     * the magnitude rather than merely acknowledging the issue. */
    printf("\n  Biomass-coupling scale check (why biomass is out of scope):\n");
    printf("  %-6s %14s %14s %16s\n", "PP", "lambda_sb", "B_l(end)", "lambda_sb*B_l");
    for (int i = 0; i < NPP; ++i) {
        double lsb = (double)lambda_sb_f(
            input.growth.lambda_sb,
            input.nutrients.nitrogen_soil_content,
            input.nutrients.phosphorus_soil_content);
        printf("  %-6d %14.8g %14.8g %16.6e\n",
            rows[i].photoperiod_h, lsb, rows[i].bl1, lsb * rows[i].bl1);
    }
    printf("\n  A stoichiometric budget would require this product to be 1.\n");

    return 0;
}

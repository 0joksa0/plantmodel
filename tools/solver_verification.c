/* ---------------------------------------------------------------------------
 * solver_verification -- numerical verification of the Ci-An coupling solver.
 *
 * Produces the evidence for the numerical-verification section of the
 * manuscript revision (reviewer issues #12, #13, #15 and the two additional
 * defects documented in docs/revision/math_audit.md):
 *
 *   results/verification/solver_residual.csv
 *       Supply-demand residual |An - gc(Ca - Ci)| for both solvers across
 *       PAR. The analytic solver satisfies the coupling to machine
 *       precision; the fixed-point scheme does not below a PAR threshold.
 *
 *   results/verification/solver_montecarlo.csv
 *       Monte Carlo census over the published sensitivity ranges
 *       (config/sensitivity_params.csv): non-convergence rate, fraction of
 *       samples returning An < 0 or Ci > 2 Ca, and the discrepancy
 *       distribution against the analytic reference.
 *
 *   results/verification/solver_parity.csv
 *       Returned An as a function of max_iter, demonstrating the period-2
 *       orbit: the value alternates with the parity of the iteration count,
 *       so the answer is an artefact of the loop bound.
 *
 *   results/verification/solver_water_stress.csv
 *       Non-convergence rate and integrated daily carbon gain as a function
 *       of leaf water potential, over a full sinusoidal 12 h photoperiod.
 *
 *   results/verification/jarvis_temperature.csv
 *       Corrected vs previous Jarvis f_T and the resulting conductance,
 *       showing that the correction is exact at the calibration temperature.
 *
 * Build: registered in CMakeLists.txt. Run from the repository root.
 * ------------------------------------------------------------------------- */
#include "model/model.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#define MC_SAMPLES 20000

/* xorshift64* -- deterministic, self-contained, so the numbers in the paper
 * are reproducible without depending on the platform's rand(). */
static unsigned long long g_rng = 20240513ULL;

static double rng_uniform(void)
{
    g_rng ^= g_rng << 13;
    g_rng ^= g_rng >> 7;
    g_rng ^= g_rng << 17;
    return (double)(g_rng >> 11) / 9007199254740992.0;
}

static double rng_range(double lo, double hi) { return lo + rng_uniform() * (hi - lo); }

static double rng_log_range(double lo, double hi)
{
    return exp(log(lo) + rng_uniform() * (log(hi) - log(lo)));
}

static int cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

static void ensure_dir(void)
{
    mkdir("results", 0755);
    mkdir("results/verification", 0755);
}

/* Baseline leaf: the default configuration of the production runs. */
static Input reference_leaf(void)
{
    Input in = generate_input();
    in.gas_exchange.leaf_temperature = REAL(25.0);
    in.gas_exchange.relative_humidity = REAL(60.0);
    in.gas_exchange.ambient_CO2_concentration = REAL(400.0);
    in.gas_exchange.leaf_water_potential = REAL(-0.6);
    in.gas_exchange.max_rubisco_carboxylation_rate = REAL(80.0);
    in.gas_exchange.respiration_rate = REAL(1.0);
    in.gas_exchange.CO2_compensation_point = REAL(42.0);
    in.gas_exchange.michaelis_constant_CO2 = REAL(404.9);
    in.gas_exchange.michaelis_constant_O2 = REAL(278.4);
    in.gas_exchange.oxygen_concentration = REAL(210.0);
    return in;
}

/* Sets gs and J for the given light level, leaving the leaf ready to solve. */
static void set_light(Input *in, double par, double jmax)
{
    double vpd = calculate_vapor_pressure_deficit(
        in->gas_exchange.leaf_temperature, in->gas_exchange.relative_humidity);
    in->gas_exchange.stomatal_conductance = calculate_stomatal_conductance_jarvis(
        par, vpd,
        in->gas_exchange.leaf_temperature,
        in->gas_exchange.ambient_CO2_concentration,
        in->gas_exchange.leaf_water_potential);
    in->gas_exchange.electron_transport_rate = compute_J(par, jmax);
}

/* --- (1) supply-demand residual across PAR ------------------------------- */
static void write_residual(void)
{
    FILE *f = fopen("results/verification/solver_residual.csv", "w");
    if (!f) { perror("solver_residual.csv"); return; }
    fprintf(f, "PAR,gs,Ci_analytic,An_analytic,residual_analytic,"
               "Ci_fixedpoint,An_fixedpoint,residual_fixedpoint,fp_converged\n");

    for (double par = 10.0; par <= 2000.0; par *= 1.15)
    {
        Input in = reference_leaf();
        set_light(&in, par, 150.0);
        double gs = in.gas_exchange.stomatal_conductance;
        double gc = gs / 1.6;
        double J = in.gas_exchange.electron_transport_rate;

        double an_a = solve_ci_analytic(&in);
        double ci_a = in.gas_exchange.intercellular_CO2;
        double res_a = fabs(an_a - gc * (in.gas_exchange.ambient_CO2_concentration - ci_a));

        in.gas_exchange.electron_transport_rate = J;
        ci_solver_reset_diagnostics();
        double an_f = iterate_ci(&in, 1000, REAL(0.0001));
        double ci_f = in.gas_exchange.intercellular_CO2;
        int conv = (ci_solver_get_diagnostics().nonconvergent == 0);
        double res_f = fabs(an_f - gc * (in.gas_exchange.ambient_CO2_concentration - ci_f));

        fprintf(f, "%.2f,%.8f,%.6f,%.6f,%.3e,%.6f,%.6f,%.3e,%d\n",
                par, gs, ci_a, an_a, res_a, ci_f, an_f, res_f, conv);
    }
    fclose(f);
    printf("wrote results/verification/solver_residual.csv\n");
}

/* --- (2) loop-parity artefact -------------------------------------------- */
static void write_parity(void)
{
    FILE *f = fopen("results/verification/solver_parity.csv", "w");
    if (!f) { perror("solver_parity.csv"); return; }
    fprintf(f, "max_iter,An_fixedpoint,Ci_fixedpoint,An_analytic\n");

    /* A water potential inside the failure band, so the map does not settle. */
    for (int m = 990; m <= 1005; m++)
    {
        Input in = reference_leaf();
        in.gas_exchange.leaf_water_potential = REAL(-1.4);
        set_light(&in, 1000.0, 150.0);
        double J = in.gas_exchange.electron_transport_rate;
        double an_f = iterate_ci(&in, m, REAL(0.0001));
        double ci_f = in.gas_exchange.intercellular_CO2;
        in.gas_exchange.electron_transport_rate = J;
        double an_a = solve_ci_analytic(&in);
        fprintf(f, "%d,%.6f,%.6f,%.6f\n", m, an_f, ci_f, an_a);
    }
    fclose(f);
    printf("wrote results/verification/solver_parity.csv\n");
}

/* --- (3) Monte Carlo census ---------------------------------------------- */
static void write_montecarlo(void)
{
    double *diff = malloc(MC_SAMPLES * sizeof(double));
    if (!diff) return;
    int nonconv = 0, negative = 0, ci_runaway = 0;
    double sum = 0.0;

    for (int i = 0; i < MC_SAMPLES; i++)
    {
        /* Ranges taken from config/sensitivity_params.csv. */
        double par  = rng_range(300.0, 2000.0);
        double temp = rng_range(10.0, 35.0);
        double rh   = rng_range(30.0, 90.0);
        double ca   = rng_range(280.0, 800.0);
        double psi  = rng_range(-1.5, -0.1);
        double vcm  = rng_log_range(20.0, 120.0);
        double jmx  = rng_log_range(50.0, 200.0);

        Input in = reference_leaf();
        in.gas_exchange.leaf_temperature = temp;
        in.gas_exchange.relative_humidity = rh;
        in.gas_exchange.ambient_CO2_concentration = ca;
        in.gas_exchange.leaf_water_potential = psi;
        in.gas_exchange.max_rubisco_carboxylation_rate = vcm;
        set_light(&in, par, jmx);
        double J = in.gas_exchange.electron_transport_rate;

        ci_solver_reset_diagnostics();
        double an_f = iterate_ci(&in, 1000, REAL(0.0001));
        double ci_f = in.gas_exchange.intercellular_CO2;
        if (ci_solver_get_diagnostics().nonconvergent) nonconv++;
        if (an_f < 0.0) negative++;
        if (ci_f > 2.0 * ca) ci_runaway++;

        in.gas_exchange.electron_transport_rate = J;
        double an_a = solve_ci_analytic(&in);

        diff[i] = fabs(an_f - an_a);
        sum += diff[i];
    }

    qsort(diff, MC_SAMPLES, sizeof(double), cmp_double);

    FILE *f = fopen("results/verification/solver_montecarlo.csv", "w");
    if (!f) { perror("solver_montecarlo.csv"); free(diff); return; }
    fprintf(f, "metric,value\n");
    fprintf(f, "samples,%d\n", MC_SAMPLES);
    fprintf(f, "nonconvergent_pct,%.3f\n", 100.0 * nonconv / MC_SAMPLES);
    fprintf(f, "negative_An_pct,%.3f\n", 100.0 * negative / MC_SAMPLES);
    fprintf(f, "Ci_exceeds_2Ca_pct,%.3f\n", 100.0 * ci_runaway / MC_SAMPLES);
    fprintf(f, "mean_abs_discrepancy,%.5f\n", sum / MC_SAMPLES);
    fprintf(f, "median_abs_discrepancy,%.5f\n", diff[MC_SAMPLES / 2]);
    fprintf(f, "p95_abs_discrepancy,%.5f\n", diff[(int)(0.95 * MC_SAMPLES)]);
    fprintf(f, "max_abs_discrepancy,%.5f\n", diff[MC_SAMPLES - 1]);
    fclose(f);
    free(diff);
    printf("wrote results/verification/solver_montecarlo.csv\n");
}

/* --- (4) water-stress dependence over a full diel cycle ------------------ */
static void write_water_stress(void)
{
    FILE *f = fopen("results/verification/solver_water_stress.csv", "w");
    if (!f) { perror("solver_water_stress.csv"); return; }
    fprintf(f, "leaf_water_potential,daylight_samples,nonconvergent_samples,"
               "nonconvergent_pct,daily_C_fixedpoint_mol_m2,"
               "daily_C_analytic_mol_m2,carbon_loss_pct\n");

    const int steps = 24 * 200;
    const double dt = 24.0 / steps;

    for (double psi = -0.1; psi >= -1.5001; psi -= 0.1)
    {
        Input in = reference_leaf();
        in.gas_exchange.leaf_water_potential = psi;
        in.photo.photoperiod = REAL(12.0);
        in.photo.max_PAR = REAL(1000.0);

        int day = 0, bad = 0;
        double c_fp = 0.0, c_an = 0.0;

        for (int i = 0; i < steps; i++)
        {
            update_light_conditions(&in, (i + 0.5) * dt);
            if (in.photo.light <= 0.0) continue;
            day++;
            set_light(&in, in.photo.light_PAR, 150.0);
            double J = in.gas_exchange.electron_transport_rate;

            ci_solver_reset_diagnostics();
            double an_f = iterate_ci(&in, 1000, REAL(0.0001));
            if (ci_solver_get_diagnostics().nonconvergent) bad++;

            in.gas_exchange.electron_transport_rate = J;
            double an_a = solve_ci_analytic(&in);

            /* Clamped as in farquhar_photosynthesis()'s RMAX(0, .). */
            c_fp += (an_f > 0.0 ? an_f : 0.0) * dt * 3600.0 * 1e-6;
            c_an += (an_a > 0.0 ? an_a : 0.0) * dt * 3600.0 * 1e-6;
        }

        double loss = (c_an > 0.0) ? 100.0 * (c_an - c_fp) / c_an : 0.0;
        fprintf(f, "%.2f,%d,%d,%.2f,%.6f,%.6f,%.2f\n",
                psi, day, bad, day ? 100.0 * bad / day : 0.0, c_fp, c_an, loss);
    }
    fclose(f);
    printf("wrote results/verification/solver_water_stress.csv\n");
}

/* --- (5) Jarvis temperature response, corrected vs previous -------------- */
static void write_jarvis_temperature(void)
{
    FILE *f = fopen("results/verification/jarvis_temperature.csv", "w");
    if (!f) { perror("jarvis_temperature.csv"); return; }
    fprintf(f, "temperature_C,f_T_corrected,f_T_previous,gs_corrected,"
               "gs_previous,gs_ratio\n");

    /* Previous implementation, reproduced verbatim from the pre-revision
     * source so the comparison is exact rather than remembered. */
    const double T0 = 25.0, T1 = 5.0, Th = 45.0;
    const double b4_old = (Th - T0) / (Th - T1);
    const double b3_old = 1.0 / ((T0 - T1) * pow(Th - T1, b4_old));
    const double b4_new = (Th - T0) / (T0 - T1);
    const double b3_new = 1.0 / ((T0 - T1) * pow(Th - T0, b4_new));
    const double gs_max = 0.70710678118654752;

    for (double T = 5.0; T <= 45.0001; T += 1.0)
    {
        double ft_new = (T > T1 && T < Th)
                            ? b3_new * (T - T1) * pow(Th - T, b4_new) : 0.0;
        if (ft_new > 1.0) ft_new = 1.0;
        double ft_old = (T > T1 && T < Th)
                            ? b3_old * (T - T1) * pow(Th - T, b4_old) : 0.0;
        if (ft_old > 1.0) ft_old = 1.0;
        if (ft_old < 0.0) ft_old = 0.0;

        double vpd = calculate_vapor_pressure_deficit(T, 60.0);
        double gs_new = calculate_stomatal_conductance_jarvis(1000.0, vpd, T, 400.0, -0.6);
        /* Undo gs_max and f_T_new, reapply f_T_old, to recover the old gs. */
        double gs_old = (ft_new > 0.0) ? gs_new / gs_max / ft_new * ft_old : 0.0;

        fprintf(f, "%.1f,%.6f,%.6f,%.8f,%.8f,%.6f\n", T, ft_new, ft_old,
                gs_new, gs_old, (gs_old > 0.0) ? gs_new / gs_old : 0.0);
    }
    fclose(f);
    printf("wrote results/verification/jarvis_temperature.csv\n");
}

int main(void)
{
    ensure_dir();
    write_residual();
    write_parity();
    write_montecarlo();
    write_water_stress();
    write_jarvis_temperature();
    return 0;
}

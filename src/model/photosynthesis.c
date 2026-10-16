#include "model/model.h"

#include <math.h>
#include <stdio.h>

static FILE *g_photo_fp = NULL;

/* ---------------------------------------------------------------------------
 * Ci-An solver state: mode selection and diagnostic counters.
 * Definitions of the accessors appear further down, next to the solvers.
 * ------------------------------------------------------------------------- */
static CiSolverDiagnostics g_ci_diag = {0, 0, 0, 0, 0};
static CiSolverMode        g_ci_mode = CI_SOLVER_ANALYTIC;

/* Fixed-point solver settings. Defaults reproduce the historical hardcoded
 * call iterate_ci(input, 1000, 1e-4) exactly: eps_an < 0 makes
 * iterate_ci_tol() apply eps_ci to both criteria. */
static real_t g_ci_eps_ci   = REAL(0.0001);
static real_t g_ci_eps_an   = REAL(-1.0);
static int    g_ci_max_iter = 1000;

void photo_dataset_open(const char *path)
{
    if (g_photo_fp)
    {
        fclose(g_photo_fp);
        g_photo_fp = NULL;
    }

    g_photo_fp = fopen(path, "w");
    if (!g_photo_fp)
    {
        perror("photo_dataset_open");
        return;
    }

    fprintf(g_photo_fp,
            "light_PAR,leaf_biomass,leaf_temperature,relative_humidity,"
            "ambient_CO2,leaf_water_potential,specific_leaf_area,"
            "VPD,stomatal_conductance,CO2_compensation_point,"
            "Vcmax,Kc,Ko,O2,J,Rd,target_An,target_photosynthesis\n");
}

void photo_dataset_close(void)
{
    if (g_photo_fp)
    {
        fflush(g_photo_fp);
        fclose(g_photo_fp);
        g_photo_fp = NULL;
    }
}

void photo_dataset_flush(void)
{
    if (g_photo_fp)
    {
        fflush(g_photo_fp);
    }
}

static void photo_dataset_log(
    Input *input,
    real_t VPD,
    real_t gs,
    real_t An,
    real_t adjusted_PH)
{
    if (!g_photo_fp)
    {
        return;
    }

    fprintf(g_photo_fp,
            "%.10Lf,%.10Lf,%.10Lf,%.10Lf,"
            "%.10Lf,%.10Lf,%.10Lf,"
            "%.10Lf,%.10Lf,%.10Lf,"
            "%.10Lf,%.10Lf,%.10Lf,%.10Lf,%.10Lf,%.10Lf,%.10Lf,%.10Lf\n",
            (long double)input->photo.light_PAR,
            (long double)input->growth.leaf_biomass,
            (long double)input->gas_exchange.leaf_temperature,
            (long double)input->gas_exchange.relative_humidity,
            (long double)input->gas_exchange.ambient_CO2_concentration,
            (long double)input->gas_exchange.leaf_water_potential,
            (long double)input->gas_exchange.specific_leaf_area,
            (long double)VPD,
            (long double)gs,
            (long double)input->gas_exchange.CO2_compensation_point,
            (long double)input->gas_exchange.max_rubisco_carboxylation_rate,
            (long double)input->gas_exchange.michaelis_constant_CO2,
            (long double)input->gas_exchange.michaelis_constant_O2,
            (long double)input->gas_exchange.oxygen_concentration,
            (long double)input->gas_exchange.electron_transport_rate,
            (long double)input->gas_exchange.respiration_rate,
            (long double)An,
            (long double)adjusted_PH);
}

real_t photosynthesis(
    real_t light,
    real_t limitation_of_photosynthetic_rate,
    real_t max_photosynthetic_rate,
    real_t nitrogen_saturation,
    real_t phosphorus_saturation,
    real_t leaf_biomass,
    real_t min_leaf_biomass)
{
    if (leaf_biomass < min_leaf_biomass)
    {
        return 0;
    }

    real_t photo = light * limitation_of_photosynthetic_rate * max_photosynthetic_rate * RMIN(nitrogen_saturation, phosphorus_saturation);
    return photo;
}

/*
 * Non-rectangular hyperbola for actual electron transport rate J from PAR.
 * Farquhar et al. (1980), Eq. 4 form widely used in practice:
 *   θ·J² - (Jmax + φ·PAR)·J + Jmax·φ·PAR = 0
 * φ = 0.385 mol e⁻/mol photons (quantum yield incl. absorptivity)
 * θ = 0.7  (curvature, Farquhar & Wong 1984)
 */
real_t compute_J(real_t PAR, real_t Jmax)
{
    const real_t phi   = REAL(0.385);
    const real_t theta = REAL(0.7);
    real_t b    = -(Jmax + phi * PAR);
    real_t disc = b * b - REAL(4.0) * theta * Jmax * phi * PAR;
    if (disc < REAL(0.0)) disc = REAL(0.0);
    return (-b - RSQRT(disc)) / (REAL(2.0) * theta);
}

real_t farquhar_photosynthesis(Input *input)
{
    if (input->photo.light <= REAL(0.0) || input->growth.leaf_biomass < input->photo.min_leaf_biomass)
    {
        return REAL(0.0);
    }
    /* The hard PAR < 20 umol m^-2 s^-1 cut-off previously applied here was a
     * discontinuity: f_PAR(20) = 0.0374 != 0, so the returned assimilation
     * jumped by roughly +0.87 umol m^-2 s^-1 across the threshold, which is a
     * hazard for any adaptive-step ODE integrator (reviewer issue #6). With
     * the analytic solver, A_n passes smoothly through zero at
     * PAR = 13.36 umol m^-2 s^-1 under default conditions, so the terminal
     * RMAX(0, .) clamp already returns exactly zero below the light
     * compensation point, continuously and without a special case. */

    input->gas_exchange.electron_transport_rate =
        compute_J(input->photo.light_PAR, input->gas_exchange.electron_transport_rate);

    real_t VPD = calculate_vapor_pressure_deficit(input->gas_exchange.leaf_temperature, input->gas_exchange.relative_humidity);
    input->gas_exchange.stomatal_conductance = calculate_stomatal_conductance_jarvis(
        input->photo.light_PAR,
        VPD,
        input->gas_exchange.leaf_temperature,
        input->gas_exchange.ambient_CO2_concentration,
        input->gas_exchange.leaf_water_potential);

    real_t An = (g_ci_mode == CI_SOLVER_FIXED_POINT)
                    ? iterate_ci_tol(input, g_ci_max_iter, g_ci_eps_ci, g_ci_eps_an)
                    : solve_ci_analytic(input);

    input->gas_exchange.convert_to_photosynthesis_per_gram_per_hour = convert_to_photosynthesis_per_gram_per_hour(
        An,
        input->growth.leaf_biomass,
        input->gas_exchange.specific_leaf_area);
    real_t adjusted_PH = input->gas_exchange.convert_to_photosynthesis_per_gram_per_hour;

    photo_dataset_log(
        input,
        VPD,
        input->gas_exchange.stomatal_conductance,
        An,
        adjusted_PH);

#ifdef PLANTMODEL_SELF_SHADING
    /* Self-shading.  convert_to_photosynthesis_per_gram_per_hour() returns a
     * per-gram rate in which every gram of leaf sees full PAR, so total
     * assimilation grows linearly with leaf mass without limit.  At day 30 the
     * 12 h plant carries 0.88 m2 of leaf at the atlas SLA, against ~0.0105 m2
     * for a real three-week Col-0 rosette (83.9x), and the residual left after
     * one shared lambda_sb is monotone in photoperiod: +21.9/-0.6/-18.9/-31.4%.
     * Both point at the same missing physics.
     *
     * Beer-Lambert over the rosette: the mean per-unit-area rate of a canopy of
     * leaf area index L is (1 - exp(-k L))/(k L), which is 1 at L -> 0 and falls
     * as 1/(kL) for a dense canopy.  L = SLA * B_leaf is dimensionless, and k is
     * the extinction coefficient -- 0.5 for a spherical leaf-angle distribution
     * (Monsi & Saeki 1953), the standard value for rosettes. */
    {
        real_t L = input->gas_exchange.specific_leaf_area * input->growth.leaf_biomass;
        real_t kL = input->photo.extinction_coefficient * L;
        if (kL > REAL(1e-9))
            adjusted_PH *= (REAL(1.0) - REXP(-kL)) / kL;
    }
#endif

#ifdef PLANTMODEL_STARCH_FEEDBACK_FARQUHAR
    /* Sink limitation.  limitation_of_photosynthetic_rate() is computed in
     * plant_compute_outputs() but only PHOTO_MODE_BASIC multiplies by it, so
     * under the Farquhar path assimilation was never down-regulated by a full
     * starch pool.  The consequence is measurable: at 12 h the pool runs to
     * 1911 against a ceiling of 72 umol C6 gFW-1, and RGR stops increasing with
     * photoperiod.  The factor is clamped at its ceiling value so it cannot
     * turn assimilation negative above the ceiling. */
    {
        real_t lim = input->photo.limitation_of_photosynthetic_rate;
        lim = RMAX(lim, input->photo.feedback_on_photosynthesis);
        lim = RMIN(lim, REAL(1.0));
        adjusted_PH *= lim;
    }
#endif

    return RMAX(REAL(0.0), adjusted_PH);
}

real_t farquhar_only_photosynthesis(Input *input)
{
    if (input->photo.light <= REAL(0.0) || input->growth.leaf_biomass < input->photo.min_leaf_biomass)
        return REAL(0.0);
    /* PAR < 20 cut-off removed; see farquhar_photosynthesis() above. */

    input->gas_exchange.electron_transport_rate =
        compute_J(input->photo.light_PAR, input->gas_exchange.electron_transport_rate);

    real_t An = (g_ci_mode == CI_SOLVER_FIXED_POINT)
                    ? iterate_ci_tol(input, g_ci_max_iter, g_ci_eps_ci, g_ci_eps_an)
                    : solve_ci_analytic(input);

    input->gas_exchange.convert_to_photosynthesis_per_gram_per_hour = convert_to_photosynthesis_per_gram_per_hour(
        An,
        input->growth.leaf_biomass,
        input->gas_exchange.specific_leaf_area);

    return RMAX(REAL(0.0), input->gas_exchange.convert_to_photosynthesis_per_gram_per_hour);
}

real_t rubisco_limited_photosynthesis(
    real_t intercellular_CO2,
    real_t CO2_compensation_point,
    real_t max_rubisco_carboxylation_rate,
    real_t michaelis_constant_CO2,
    real_t michaelis_constant_O2,
    real_t oxygen_concentration)
{
    real_t denominator = intercellular_CO2 + (michaelis_constant_CO2 * (REAL(1.0) + (oxygen_concentration / michaelis_constant_O2)));
    return max_rubisco_carboxylation_rate * (intercellular_CO2 - CO2_compensation_point) / denominator;
}

real_t light_limited_photosynthesis(
    real_t intercellular_CO2,
    real_t CO2_compensation_point,
    real_t electron_transport_rate)
{
    real_t denominator = (REAL(4.5) * intercellular_CO2) + (REAL(10.5) * CO2_compensation_point);
    return electron_transport_rate * (intercellular_CO2 - CO2_compensation_point) / denominator;
}

real_t net_photosynthesis(
    real_t rubisco_limited_rate,
    real_t light_limited_rate,
    real_t respiration_rate)
{
    real_t limiting_rate = RMIN(rubisco_limited_rate, light_limited_rate);
    return limiting_rate - respiration_rate;
}

real_t convert_to_photosynthesis_per_gram_per_hour(
    real_t net_photosynthesis_rate,
    real_t leaf_biomass,
    real_t specific_leaf_area)
{
    /* PH [umol C6 gFW-1 h-1] = An [umol CO2 m-2 s-1] * SLA [m2 gFW-1] * 3600/6.
     * No free constant: at the atlas SLA 0.005500 m2 gFW-1 (Quantitative Plant
     * Biology 10.1017/qpb.2024.1, Table 1) this gives 3.3000, against the
     * independently measured requirement 3.3015.  base_scale = 1.10 and
     * ref_specific_leaf_area = 0.025 were both unsourced and are removed. */
    real_t safe_sla = RMAX(specific_leaf_area, REAL(1e-6));

    (void)leaf_biomass;
    return net_photosynthesis_rate * safe_sla * REAL(600.0);
}

real_t intercellular_CO2(
    real_t ambient_CO2_concentration,
    real_t net_photosynthesis_rate,
    real_t stomatal_conductance)
{
    return ambient_CO2_concentration - (net_photosynthesis_rate / stomatal_conductance);
}

void ci_solver_reset_diagnostics(void)
{
    g_ci_diag.calls = 0;
    g_ci_diag.nonconvergent = 0;
    g_ci_diag.clamped_ci = 0;
    g_ci_diag.closed_stomata = 0;
    g_ci_diag.iterations_total = 0;
}

CiSolverDiagnostics ci_solver_get_diagnostics(void) { return g_ci_diag; }

void ci_solver_set_mode(CiSolverMode mode) { g_ci_mode = mode; }
CiSolverMode ci_solver_get_mode(void) { return g_ci_mode; }

void ci_solver_set_tolerances(real_t eps_ci, real_t eps_an)
{
    g_ci_eps_ci = eps_ci;
    g_ci_eps_an = eps_an;
}

void ci_solver_get_tolerances(real_t *eps_ci, real_t *eps_an)
{
    if (eps_ci) *eps_ci = g_ci_eps_ci;
    if (eps_an) *eps_an = g_ci_eps_an;
}

int  ci_solver_get_max_iter(void)          { return g_ci_max_iter; }
void ci_solver_set_max_iter(int max_iter)  { g_ci_max_iter = (max_iter > 0) ? max_iter : 1; }

/* ---------------------------------------------------------------------------
 * solve_ci_analytic(): closed-form solution of the supply-demand coupling.
 *
 * Supply (Fick):   An = gc * (Ca - Ci),          gc = gs / 1.6
 * Demand (FvCB):   An = num * (Ci - Gstar) / (d1*Ci + d0) - Rd
 *
 * Equating and clearing the denominator,
 *   gc*(Ca - Ci)*(d1*Ci + d0) = num*(Ci - Gstar) - Rd*(d1*Ci + d0),
 * gives a*Ci^2 + b*Ci + c = 0 with
 *   a = -gc*d1
 *   b =  gc*(d1*Ca - d0) - num + Rd*d1
 *   c =  gc*d0*Ca + num*Gstar + Rd*d0
 *
 * Solved for each limitation branch; the admissible root is substituted back
 * into the demand function and the co-limiting (minimum) branch is returned,
 * matching net_photosynthesis()'s RMIN(Ac, Aj) - Rd.
 *
 * Root admissibility is Ci >= 0 only. There is deliberately NO upper bound.
 * Requiring Ci <= Ca is wrong below the light compensation point: there An < 0
 * because respiration exceeds gross assimilation, so Fick's law An = gc(Ca-Ci)
 * demands Ci > Ca -- CO2 diffuses OUT of the leaf. An earlier version of this
 * function rejected roots above Ca*1.001 and fell through to the "no admissible
 * root" branch, which returned Ci = 0 and hence a spuriously POSITIVE An from
 * the wrong branch of the demand curve. Under production settings that made
 * A_n wrong for all PAR below about 27 umol m^-2 s^-1 -- a regime every diel
 * run passes through twice a day, at dawn and at dusk.
 *
 * Exactly one root per branch is admissible, and no search is needed to know
 * which: a = -gc*d1 < 0 and c = gc*d0*Ca + num*Gstar + Rd*d0 > 0, so the
 * product of the roots c/a is negative and each branch has precisely one
 * positive root. Selecting the minimum An over the two branches therefore
 * picks the true solution of the coupled system, not merely a candidate.
 *
 * The gc -> 0 limit is handled physically: closed stomata admit no CO2, so
 * assimilation cannot proceed and the leaf respires, An = -Rd. The previous
 * implementation clamped gs at 1e-6 and returned An = 0 there, which both
 * discarded night/low-light respiration and drove the period-2 orbit.
 * ------------------------------------------------------------------------- */
real_t solve_ci_analytic(Input *input)
{
    const real_t Ca    = input->gas_exchange.ambient_CO2_concentration;
    const real_t Gstar = input->gas_exchange.CO2_compensation_point;
    const real_t Rd    = input->gas_exchange.respiration_rate;
    const real_t Vcmax = input->gas_exchange.max_rubisco_carboxylation_rate;
    const real_t Jrate = input->gas_exchange.electron_transport_rate;
    const real_t Kc    = input->gas_exchange.michaelis_constant_CO2;
    const real_t Ko    = input->gas_exchange.michaelis_constant_O2;
    const real_t O2    = input->gas_exchange.oxygen_concentration;

    g_ci_diag.calls++;

    const real_t gc = input->gas_exchange.stomatal_conductance / REAL(1.6);
    if (gc <= REAL(0.0))
    {
        g_ci_diag.closed_stomata++;
        input->gas_exchange.intercellular_CO2     = Gstar;
        input->gas_exchange.net_photosynthesis    = -Rd;
        input->gas_exchange.net_photosynthesis_rate = -Rd;
        return -Rd;
    }

    /* Branch coefficients: {numerator, d1, d0} for Ac and Aj respectively. */
    const real_t branch[2][3] = {
        { Vcmax, REAL(1.0), Kc * (REAL(1.0) + O2 / Ko) },
        { Jrate, REAL(4.5), REAL(10.5) * Gstar }
    };

    real_t best_An = REAL(0.0);
    real_t best_Ci = REAL(0.0);
    int    found   = 0;

    for (int k = 0; k < 2; k++)
    {
        const real_t num = branch[k][0];
        const real_t d1  = branch[k][1];
        const real_t d0  = branch[k][2];

        const real_t a = -gc * d1;
        const real_t b = gc * (d1 * Ca - d0) - num + Rd * d1;
        const real_t c = gc * d0 * Ca + num * Gstar + Rd * d0;

        real_t disc = b * b - REAL(4.0) * a * c;
        if (disc < REAL(0.0)) continue;
        disc = RSQRT(disc);

        const real_t roots[2] = { (-b + disc) / (REAL(2.0) * a),
                                  (-b - disc) / (REAL(2.0) * a) };

        for (int r = 0; r < 2; r++)
        {
            const real_t Ci = roots[r];
            /* Ci >= 0 is the only constraint; Ci > Ca is physical whenever
             * An < 0 (net efflux below the light compensation point). */
            if (Ci < REAL(0.0)) continue;

            const real_t An = num * (Ci - Gstar) / (d1 * Ci + d0) - Rd;
            if (!found || An < best_An)
            {
                best_An = An;
                best_Ci = Ci;
                found   = 1;
            }
        }
    }

    if (!found)
    {
        /* No admissible root: supply cannot meet demand anywhere in [0, Ca].
         * The physical limit is a fully drawn-down intercellular space. */
        g_ci_diag.clamped_ci++;
        best_Ci = REAL(0.0);
        best_An = -Rd;
    }

    input->gas_exchange.intercellular_CO2       = best_Ci;
    input->gas_exchange.net_photosynthesis      = best_An;
    input->gas_exchange.net_photosynthesis_rate = best_An;
    return best_An;
}

real_t iterate_ci(Input *input, int max_iter, real_t epsilon)
{
    return iterate_ci_tol(input, max_iter, epsilon, REAL(-1.0));
}

/* Legacy fixed-point scheme. Retained for the convergence study only.
 * eps_ci applies to |dCi| (units umol mol^-1), eps_an to |dAn| (units
 * umol m^-2 s^-1). Passing eps_an <= 0 applies eps_ci to both, which
 * reproduces the original single-tolerance behaviour exactly. */
real_t iterate_ci_tol(Input *input, int max_iter, real_t eps_ci, real_t eps_an)
{
    real_t Ci = input->gas_exchange.ambient_CO2_concentration * REAL(0.7);
    real_t prev_An = 0.0;
    real_t An = 0.0;
    real_t gs = RMAX(input->gas_exchange.stomatal_conductance, REAL(1e-6));
    int    converged = 0;
    int    used_iter = 0;

    if (eps_an <= REAL(0.0)) eps_an = eps_ci;
    g_ci_diag.calls++;

    for (int i = 0; i < max_iter; i++)
    {
        used_iter = i + 1;
        real_t Ac = rubisco_limited_photosynthesis(
            Ci,
            input->gas_exchange.CO2_compensation_point,
            input->gas_exchange.max_rubisco_carboxylation_rate,
            input->gas_exchange.michaelis_constant_CO2,
            input->gas_exchange.michaelis_constant_O2,
            input->gas_exchange.oxygen_concentration);

        real_t Aj = light_limited_photosynthesis(
            Ci,
            input->gas_exchange.CO2_compensation_point,
            input->gas_exchange.electron_transport_rate);

        An = net_photosynthesis(Ac, Aj, input->gas_exchange.respiration_rate);

        /* gc = gs/1.6: CO2 diffusivity is ~1.6x lower than H2O in air */
        real_t gc = gs / REAL(1.6);
        real_t new_Ci = input->gas_exchange.ambient_CO2_concentration - An / gc;
        new_Ci = RMAX(REAL(0.0), new_Ci);

        if (RABS(new_Ci - Ci) < eps_ci && RABS(An - prev_An) < eps_an)
        {
            converged = 1;
            Ci = new_Ci;
            break;
        }

        Ci = new_Ci;
        prev_An = An;
    }

    g_ci_diag.iterations_total += used_iter;
    if (!converged) g_ci_diag.nonconvergent++;

    input->gas_exchange.intercellular_CO2 = Ci;
    input->gas_exchange.net_photosynthesis = An;
    input->gas_exchange.net_photosynthesis_rate = An;
    return An;
}

real_t calculate_stomatal_conductance_jarvis(
    real_t incoming_PAR,
    real_t vapor_pressure_deficit,
    real_t leaf_temperature_celsius,
    real_t ambient_CO2_concentration,
    real_t leaf_water_potential)
{
    const real_t b1 = REAL(0.6);
    const real_t b2 = REAL(0.002);
    const real_t b10 = REAL(0.05);
    real_t q = b10 / b1;
    real_t f_PAR = (incoming_PAR <= q) ? b10 : (b1 * b2 * (incoming_PAR - q)) / (b1 + b2 * (incoming_PAR - q));

    const real_t b5 = REAL(0.2);
    real_t f_VPD = RMAX(REAL(0.0), REAL(1.0) - b5 * vapor_pressure_deficit);

    /* Jarvis (1976) / Stewart (1988) temperature response.
     *
     * f_T = b3 * (T - T1) * (Th - T)^b4, normalised so that f_T(T0) = 1 with
     * the maximum located at T0 and zeros at T1 and Th.
     *
     * Two conditions fix b3 and b4:
     *   stationarity at T0  =>  b4 = (Th - T0) / (T0 - T1)
     *   f_T(T0) = 1         =>  b3 = 1 / ((T0 - T1) * (Th - T0)^b4)
     *
     * The previous implementation used b4 = (Th - T0)/(Th - T1) = 0.5 and
     * (Th - T1)^b4 in the normaliser, which gave f_T(25) = 0.7071 with the
     * peak displaced to 31.7 C. Both are corrected here; see
     * docs/revision/math_audit.md section 4.2. The RMIN(.., 1.0) clamp is no
     * longer needed because the normalised form cannot exceed unity. */
    const real_t T0 = REAL(25.0);
    const real_t T1 = REAL(5.0);
    const real_t Th = REAL(45.0);
    const real_t b4 = (Th - T0) / (T0 - T1);
    const real_t b3 = REAL(1.0) / ((T0 - T1) * RPOW((Th - T0), b4));
    real_t f_T = REAL(0.0);
    if (leaf_temperature_celsius > T1 && leaf_temperature_celsius < Th)
    {
        f_T = b3 * (leaf_temperature_celsius - T1) * RPOW((Th - leaf_temperature_celsius), b4);
    }
    f_T = RMAX(REAL(0.0), RMIN(f_T, REAL(1.0)));

    const real_t b8 = REAL(0.2);
    real_t f_CO2 = REAL(1.0);
    if (ambient_CO2_concentration < 100.0)
    {
        f_CO2 = 1.0;
    }
    else if (ambient_CO2_concentration < 1000.0)
    {
        const real_t b7 = (1.0 - b8) / 900.0;
        f_CO2 = 1.0 - b7 * (ambient_CO2_concentration - 100.0);
    }
    else
    {
        f_CO2 = b8;
    }

    const real_t b6 = 3.0;
    const real_t psi_min = -1.5;
    real_t sigma_psi = leaf_water_potential - psi_min;
    real_t f_water = 1.0 - REXP(-b6 * sigma_psi);
    f_water = RMAX(0.0, RMIN(f_water, 1.0));

    /* Maximum (unstressed) stomatal conductance.
     *
     * The conventional Jarvis form is gs = gs_max * PROD f_i with every
     * f_i in [0,1] dimensionless. The previous implementation had no gs_max
     * term at all: its role was absorbed by b1 = 0.6 inside f_PAR, which
     * simultaneously sets the PAR half-saturation shape, conflating two
     * physically distinct roles in one fitted constant.
     *
     * The value below is exactly the mis-normalised f_T(25 C) of the previous
     * implementation, 1/sqrt(2). Because all production runs hold
     * T = 25 C fixed (config/default.config), factoring it out here makes the
     * f_T correction above an exact reparametrisation: conductance at the
     * calibration temperature is reproduced to machine precision, so no
     * recalibration is required for this change. Only the response away from
     * 25 C differs, which is the intended correction. See
     * docs/revision/math_audit.md section 4.2. */
    const real_t gs_max = REAL(0.70710678118654752);

    return gs_max * f_PAR * f_VPD * f_T * f_CO2 * f_water;
}

real_t calculate_vapor_pressure_deficit(real_t temperature_celsius, real_t relative_humidity_percent)
{
    real_t saturation_vapor_pressure = REAL(0.61078) * REXP((REAL(17.27) * temperature_celsius) / (temperature_celsius + REAL(237.3)));
    real_t actual_vapor_pressure = saturation_vapor_pressure * (relative_humidity_percent / REAL(100.0));
    return saturation_vapor_pressure - actual_vapor_pressure;
}

void update_light_conditions(Input *input, real_t current_hour)
{
    real_t sunrise = 0.0;
    real_t sunset = sunrise + input->photo.photoperiod;

    int is_day = (current_hour >= sunrise && current_hour < sunset);

    if (is_day)
    {
        input->photo.light = 1.0;
        if (input->photo.light_profile > REAL(0.5))
        {
            /* Square wave: constant irradiance across the light period, as a
             * growth chamber delivers it.  The daily light integral is then
             * max_PAR * photoperiod rather than (2/pi) * max_PAR * photoperiod,
             * so a run intended to match a chamber experiment must set max_PAR
             * to the chamber irradiance, not to a peak value. */
            input->photo.light_PAR = input->photo.max_PAR;
        }
        else
        {
            real_t day_fraction = (current_hour - sunrise) / (sunset - sunrise);
            input->photo.light_PAR = input->photo.max_PAR * RSIN(day_fraction * M_PI);
        }
    }
    else
    {
        input->photo.light = 0.0;
        input->photo.light_PAR = 0.0;
    }

    /* Thermoperiod, applied only when the caller has supplied both endpoints. */
    if (input->gas_exchange.leaf_temperature_day   > TEMPERATURE_UNSET &&
        input->gas_exchange.leaf_temperature_night > TEMPERATURE_UNSET)
    {
        input->gas_exchange.leaf_temperature = is_day
            ? input->gas_exchange.leaf_temperature_day
            : input->gas_exchange.leaf_temperature_night;
    }
}

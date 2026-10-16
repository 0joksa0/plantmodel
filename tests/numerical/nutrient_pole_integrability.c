/*
 * Regression: the adaptive step must not collapse on the nutrient pole.
 *
 * History.  Two parameter regions were reported as divergent and were excluded
 * from every optimiser and sweep:
 *
 *   (1) 12 h photoperiod, chamber configuration, kappa = 3.25, lambda_sdr < 2;
 *   (2) 12 h photoperiod, production configuration, lambda_sdi near 0.5163.
 *
 * Neither was a divergence.  In both, the nutrient pools N and P are driven to
 * zero, and the imbalance ratio (N - R*P)/(N + R*P + eps) shared by
 * nitrogen_affinity() and phosphorus_affinity() -- together with
 * stoichiometric_signal() -- has a pole on the ray N = -R*P, which passes
 * through the origin of the admissible quadrant.  A Runge-Kutta trial stage
 * that takes a near-empty pool slightly negative lands on it, the embedded
 * error estimate cannot be met at any step size, and the controller drives dt
 * toward zero (8.6e-157 observed) while t stands still.  Every state stays
 * finite, so the NaN guard never fires; the step budget eventually trips and
 * the run is labelled divergent.
 *
 * The fix takes the positive part of each pool inside the saturating ratios
 * only (nutrient_ratio_pool() in src/model/nutrients.c).  This is the identity
 * wherever the model is defined, so the tests below assert two things at once:
 * that the previously stalling cases now integrate to t_end quickly, and that
 * a case which always worked is unchanged.
 *
 * What this test does NOT assert: that the trajectories are biologically
 * sensible.  They are not -- N reaches about -30 in case (1), i.e. the state
 * leaves the domain the model is valid on.  That is now reported separately
 * through result.left_domain, and is checked here as a property of the
 * diagnosis rather than being silently tolerated.
 */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "config/config.h"
#include "model/model.h"
#include "model/state.h"
#include "simulation/simulation.h"

int gui_main_thread(void* a) { (void)a; return 0; }

static int failures = 0;

static void check(int ok, const char* what)
{
    printf("  [%s] %s\n", ok ? "PASS" : "FAIL", what);
    if (!ok)
        failures++;
}

typedef struct {
    int    diverged;
    int    stalled;
    int    left_domain;
    double domain_min;
    int    filled;
    int    total;
    double rgr;
    double seconds;
} Outcome;

/* photoperiod-specific settings used by the recalibration study; kappa is
 * applied through the SLA, exactly as convert_to_photosynthesis_per_gram_per_hour
 * scales by base_scale * SLA / 0.025. */
static Outcome run_case(double photoperiod, double lambda_sdr, double lambda_sdi,
                        double kappa, int chamber, int days)
{
    SimulationConfig cfg = simulation_default_config();
    Input in = generate_input();
    if (config_load_file("config/default.config", &cfg, &in) != 0) {
        fprintf(stderr, "cannot load config/default.config (run from the repo root)\n");
        Outcome bad = { 0 }; bad.diverged = -1; return bad;
    }
    cfg.days = days;
    cfg.gui_enabled = 0;
    cfg.write_csv = 0;

    in.photo.photo_mode = (int)PHOTO_MODE_FARQUHAR_JARVIS;
    in.photo.photoperiod = (real_t)photoperiod;
    in.carbohydrates.lambda_sdr = (real_t)lambda_sdr;
    in.carbohydrates.lambda_sdi = (real_t)lambda_sdi;
    in.gas_exchange.specific_leaf_area = (real_t)(0.025 * (kappa / 1.10));

    if (chamber) {
        in.photo.light_profile = REAL(1.0);
        in.photo.max_PAR = REAL(160.0);
        in.gas_exchange.leaf_temperature_day = REAL(21.0);
        in.gas_exchange.leaf_temperature_night = REAL(19.0);
        in.gas_exchange.leaf_temperature = REAL(21.0);
    }
    if (photoperiod == 12.0) {
        in.photo.feedback_on_photosynthesis = REAL(0.62);
        in.carbohydrates.lambda_sni = REAL(0.08);
        in.growth.lambda_sb = REAL(0.00585);
    }

    const real_t fw0 = in.growth.total_biomass;
    SimulationResult r;
    memset(&r, 0, sizeof(r));

    const clock_t t0 = clock();
    simulate_days(&cfg, &in, &r);
    const double secs = (double)(clock() - t0) / (double)CLOCKS_PER_SEC;

    Outcome o;
    o.diverged    = r.diverged;
    o.stalled     = r.stalled;
    o.left_domain = r.left_domain;
    o.domain_min  = r.domain_min;
    o.filled      = r.filled_steps;
    o.total       = r.total_steps;
    o.rgr         = (double)compute_rgr(fw0, in.growth.total_biomass, (real_t)(days - 1));
    o.seconds     = secs;
    simulation_result_free(&r);
    return o;
}

int main(void)
{
    /* Wall-clock ceiling.  Pre-fix these cases burned the full 10M-step budget
     * (12-15 s each on the reference machine); post-fix they complete in well
     * under a second.  One second is loose enough not to be flaky on a loaded
     * machine and tight enough that a returning stall fails the test rather
     * than merely slowing it down. */
    const double MAX_SECONDS = 1.0;

    printf("Case 1: 12 h, chamber, kappa = 3.25, lambda_sdr = 1.0 (was: step budget exhausted)\n");
    {
        Outcome o = run_case(12.0, 1.0, 0.10, 3.25, 1, 30);
        printf("    diverged=%d stalled=%d left_domain=%d min_pool=%.4g filled=%d/%d rgr=%.6g %.3f s\n",
               o.diverged, o.stalled, o.left_domain, o.domain_min, o.filled, o.total, o.rgr, o.seconds);
        check(o.diverged == 0, "integrates to t_end without exhausting the step budget");
        check(o.stalled == 0, "not classified as stalled");
        check(o.filled == o.total, "every sample slot filled");
        check(isfinite(o.rgr) && o.rgr > 0.0, "finite positive RGR");
        check(o.seconds < MAX_SECONDS, "completes well inside the pre-fix stall time");
        /* The trajectory is integrable but not admissible: N goes negative.
         * Asserting this keeps the distinction honest -- if a later change
         * makes the pools stay positive, this line should be revisited, not
         * quietly deleted. */
        check(o.left_domain == 1, "reports that a nutrient pool left the admissible quadrant");
    }

    printf("Case 2: 12 h, production, lambda_sdi at the bracketed threshold (was: step budget exhausted)\n");
    for (int k = 0; k < 2; ++k) {
        const double lsdi = (k == 0) ? 0.5162 : 0.5165;
        Outcome o = run_case(12.0, 0.25, lsdi, 1.10, 0, 8);
        printf("    lambda_sdi=%.4f diverged=%d stalled=%d left_domain=%d min_pool=%.4g filled=%d/%d %.3f s\n",
               lsdi, o.diverged, o.stalled, o.left_domain, o.domain_min, o.filled, o.total, o.seconds);
        check(o.diverged == 0, "integrates to t_end");
        check(o.filled == o.total, "every sample slot filled");
        check(o.seconds < MAX_SECONDS, "completes well inside the pre-fix stall time");
    }

    printf("Case 3: control -- a configuration that always integrated must be unaffected\n");
    {
        Outcome o = run_case(12.0, 2.0, 0.10, 3.25, 1, 30);
        printf("    diverged=%d stalled=%d left_domain=%d filled=%d/%d rgr=%.12g %.3f s\n",
               o.diverged, o.stalled, o.left_domain, o.filled, o.total, o.rgr, o.seconds);
        check(o.diverged == 0, "still integrates");
        check(o.left_domain == 0, "pools stay in the admissible quadrant");
        /* Pinned from a pre-fix binary built with the project's own flags
         * (-O3 -ffast-math, SOLVER_REAL_IS_DOUBLE).  The fix is the identity
         * in this region, so any drift means the regularisation has escaped
         * it.  Note the value is flag-sensitive: without
         * -DSOLVER_REAL_IS_DOUBLE the RMAX macro resolves to the float
         * overload and this reads 0.285223895074 instead. */
        check(fabs(o.rgr - 0.28522390262) < 1e-9, "RGR bit-stable against the pre-fix value");
    }

    printf("\n%s (%d failure%s)\n", failures ? "FAILED" : "OK", failures, failures == 1 ? "" : "s");
    return failures ? 1 : 0;
}
#include "model/model.h"

/* Starch feedback on assimilation.
 *
 * The law is linear in starch and, as written, unbounded below: it returns
 * feedback_on_photosynthesis exactly at starch = max_starch, which is the
 * intent, but continues along the same line above it and crosses zero near
 * starch/max_starch = 1/(1-f).  A negative limitation reverses the sign of
 * assimilation, which is not a physiological state; the adaptive step then
 * collapses while every state stays finite, so the non-finite guard never
 * fires.  Simulated starch does exceed the ceiling (1.30x at 6 h, 1.48x at
 * 8 h, 3.01x at 12 h), so this branch is reached in ordinary runs.
 *
 * PLANTMODEL_STARCH_FEEDBACK_BOUNDED clamps the return to its value at the
 * ceiling.  The clamp is the identity wherever starch <= max_starch, i.e.
 * over the whole range in which the law was derived, so it cannot change a
 * result that was valid to begin with. */
real_t limitation_of_photosynthetic_rate(
    real_t starch,
    real_t feedback_on_photosynthesis,
    real_t max_starch)
{
    real_t lim = feedback_on_photosynthesis + ((REAL(1.0) - feedback_on_photosynthesis) * ((max_starch - starch) / max_starch));
#ifdef PLANTMODEL_STARCH_FEEDBACK_BOUNDED
    lim = RMAX(lim, feedback_on_photosynthesis);
#endif
    return lim;
}

real_t max_starch(
    real_t max_starch_degradation_rate,
    real_t photoperiod)
{

    return max_starch_degradation_rate * (REAL(24.0) - photoperiod);
}

real_t starch_degradation(
    real_t max_starch_degradation_rate,
    real_t max_sucrose,
    real_t starch_night_start,
    real_t min_starch,
    real_t sucrose,
    real_t light,
    real_t time,
    real_t photoperiod)
{
    if (starch_night_start < min_starch) {
        return 0;
    }
    real_t degradation_rate = ((time / REAL(24.0)) + (REAL(1.0) - (time / REAL(24.0))) * (REAL(1.0) - (sucrose / (sucrose + max_sucrose)))) * ((starch_night_start - min_starch) / (REAL(24.0) - photoperiod));

    return (REAL(1.0) - light) * RMIN(max_starch_degradation_rate, degradation_rate);
}

real_t starch_sucrose_partition(
    real_t light,
    real_t starch_partition_coeff,
    real_t min_sucrose,
    real_t max_sucrose,
    real_t sucrose,
    real_t lambda_sdr,
    real_t lambda_sdi,
    real_t lambda_sni)
{
    real_t s_min_fract = (min_sucrose / (min_sucrose + sucrose));
    real_t s_max_fract = (sucrose / (sucrose + max_sucrose));
    real_t day_part = light * ((-starch_partition_coeff * lambda_sdr * s_min_fract) + ((REAL(1.0) - starch_partition_coeff) * lambda_sdi * s_max_fract));

    real_t night_part = (REAL(1.0) - light) * (REAL(1.0) - starch_partition_coeff) * lambda_sni * s_min_fract;

    return day_part + night_part;
}

real_t starch_sucrose_partition_f(
    real_t t,
    real_t starch_partition_coeff,
    void* params)
{
    Input* input = (Input*)params;
    (void)t;

    return starch_sucrose_partition(
        input->photo.light,
        starch_partition_coeff,
        input->carbohydrates.min_sucrose,
        input->carbohydrates.max_sucrose,
        input->carbohydrates.sucrose,
        input->carbohydrates.lambda_sdr,
        input->carbohydrates.lambda_sdi,
        input->carbohydrates.lambda_sni);
}

real_t starch_production(
    real_t photosynthesis,
    real_t starch_partition_coeff,
    real_t starch_degradation_rate)
{

    return (starch_partition_coeff * photosynthesis) - starch_degradation_rate;
}

real_t starch_production_f(
    real_t t,
    real_t starch,
    void* params)
{
    Input* input = (Input*)params;
    (void)t;
    (void)starch;
    return starch_production(input->photo.photosynthesis, input->carbohydrates.starch_partition_coeff, input->carbohydrates.starch_degradation_rate);
}

real_t sucrose_production(
    real_t sucrose_partition_coeff,
    real_t photosynthesis,
    real_t starch_degradation_rate,
    real_t uptake_cost,
    real_t transport_cost,
    real_t respiration_frequency,
    real_t sucrose,
    real_t sucrose_loading_frequency,
    real_t night_efficiency_starch,
    real_t light)
{
    return ((REAL(1.0) - sucrose_partition_coeff) * photosynthesis) + starch_degradation_rate - uptake_cost - transport_cost - maintenance_respiration_share(respiration_frequency, sucrose, light) - (sucrose_loading_frequency * night_efficiency_starch);
}

real_t sucrose_production_f(
    real_t t,
    real_t sucrose,
    void* params)
{
    Input* in = (Input*)params;
    (void)t;

    return sucrose_production(in->carbohydrates.starch_partition_coeff,
        in->photo.photosynthesis,
        in->carbohydrates.starch_degradation_rate,
        in->carbohydrates.uptake_cost,
        in->carbohydrates.transport_cost,
        in->carbohydrates.respiration_frequency,
        sucrose,
        in->carbohydrates.sucrose_loading_frequency,
        in->carbohydrates.night_efficiency_starch,
        in->photo.light);
}

real_t night_efficiency_starch(
    real_t sucrose,
    real_t max_sucrose,
    real_t lambda_g,
    real_t light,
    real_t starch_partition_coeff,
    real_t sink_activation_half,
    real_t sink_activation_hill)
{
#ifndef PLANTMODEL_SINK_ACTIVATION
    (void)sink_activation_half;
    (void)sink_activation_hill;
#endif

    real_t available = RMIN(sucrose, max_sucrose);
#ifdef PLANTMODEL_SINK_ACTIVATION
    /* Sugar-induced sink activation.
     *
     * As written, export and growth are both phi_load * nes with nes linear in
     * sucrose, so the drain is first order and the steady-state pool is
     * source/k_eff: the pool must track the source.  Measured, it does not.
     * Whole-tissue sucrose is flat across the photoperiod series -- 2.397,
     * 2.293, 2.774, 2.523 umol gFW-1, a spread of 1.05x -- while daily
     * assimilation rises 2.24x.  The model spreads 2.61x over the same series.
     * Sucrose is held near a set point, which a constant-coefficient first-order
     * drain cannot do.
     *
     * Sucrose induces its own consumption: it raises transcription of cell-wall
     * and vacuolar invertases and of sink-tissue transporters, so utilisation
     * capacity rises with the pool rather than staying proportional to it.  The
     * standard representation is a cooperative (Hill) activation factor on the
     * sink term, giving an effective export exponent above one and therefore a
     * pool that grows sublinearly with supply.
     *
     * activation = s^h / (K^h + s^h), h the Hill coefficient and K the
     * half-activation sucrose concentration.  At s << K the factor vanishes and
     * the pool fills; at s >> K it approaches one and previous behaviour
     * returns. */
    {
        /* The saturating Hill form s^h/(K^h+s^h) is bounded by one, so it can
         * only WEAKEN the drain and the pool grows further -- measured, it moved
         * unscaled sucrose nRMSE from 0.805 to 0.810.  Sink activation has to
         * strengthen the drain as the pool fills, so the multiplier must exceed
         * one above the set point: 1 + (s/K)^h.  The drain then scales as
         * s^(1+h) at high s, the steady-state pool grows as source^(1/(1+h)),
         * and sucrose is held near K while supply varies. */
        real_t r = sucrose / sink_activation_half;
        available *= REAL(1.0) + RPOW(r, sink_activation_hill);
    }
#endif
    return available * (lambda_g + ((REAL(1.0) - lambda_g) * (REAL(1.0) - light + (light * (REAL(1.0) - starch_partition_coeff)))));
}

real_t uptake_cost(
    real_t nitrogen_uptake_sucrose_consumption,
    real_t nitrogen_uptake,
    real_t phosphorus_uptake_sucrose_consumption,
    real_t phosphorus_uptake)
{
    return (nitrogen_uptake_sucrose_consumption * nitrogen_uptake) + (phosphorus_uptake_sucrose_consumption * phosphorus_uptake);
}

/* ---------------------------------------------------------------------------
 * growth_cost(): the carbon price of building biomass.
 *
 * WHY THIS TERM HAS TO EXIST.  Integrating the two carbohydrate ODEs over a
 * periodic diel cycle gives, exactly,
 *
 *     INT A dt = INT [ uc + tc + r_m*S + lambda_L*nes ] dt
 *
 * i.e. the four sucrose sinks already account for ALL gross assimilation, while
 *
 *     leaf_growth = lambda_sb * (1-alpha) * lambda_L * nes * B_leaf
 *
 * credits biomass from the same lambda_L*nes flux with an EXTRA factor of
 * biomass.  Growth is driven by the phloem-load flux but never charged against
 * it, so the carbon budget closes at 1.000 no matter how much mass is built.
 * That is why raising the conversion factor makes the system stiff: the extra
 * carbon has nowhere to go except the pools.
 *
 * THE STOICHIOMETRY IS MEASURED, NOT FITTED.  Building one gram of dry tissue
 * consumes Y_G grams of glucose beyond the carbon incorporated (growth
 * respiration); Penning de Vries et al. give Y_G = 1.4 g glucose gDW-1 for leaf
 * tissue, and Sulpice et al. (2014) measure protein synthesis alone at 21-30%
 * of daily carbon gain -- a term the model had no representation for.
 *
 * Converting to the model's units, with M(hexose) = 180.156 g mol-1 and the
 * FW:DW ratio phi,
 *
 *     dB/dt [gFW h-1]  costs  dB/dt * Y_G / (M * 1e-6 * phi)  [umol C6 gFW-1 h-1]
 *
 * so growth_stoichiometry = Y_G / (M*1e-6*phi) is the umol of hexose charged
 * per gram of fresh biomass produced.  It is a sourced constant, not a knob. */
/* Maintenance respiration of the shoot.
 *
 * The Farquhar solver returns NET assimilation: day respiration R_d is already
 * subtracted inside solve_ci_analytic().  Charging r_m * S over the whole
 * sucrose pool in every hour therefore bills leaf respiration a second time
 * during the light period.  Measured respiration is 9.8-14.9% of daily carbon
 * (Sulpice et al. 2014) and 13.8% from the leaf-respiration series in
 * data/paper/N_data.csv, against 23-26% in this model.
 *
 * The term is therefore charged only when the lights are off, where no
 * net-assimilation credit has been taken.  Measured against the model's own
 * carbon budget this moves the respiratory share from 22.7-25.7% to
 * 9.7-17.1% of daily assimilation, i.e. from three to five times the measured
 * value into the measured range.  PLANTMODEL_RESPIRATION_ALL_HOURS restores
 * the previous behaviour for the ablation study. */
real_t maintenance_respiration_share(real_t respiration_frequency,
                                     real_t sucrose,
                                     real_t light)
{
#ifdef PLANTMODEL_RESPIRATION_ALL_HOURS
    (void)light;
    return respiration_frequency * sucrose;      /* double-charges R_d in the light */
#else
    return (REAL(1.0) - light) * respiration_frequency * sucrose;
#endif
}

real_t growth_cost(
    real_t growth_stoichiometry,
    real_t leaf_growth_rate,
    real_t root_growth_rate)
{
    real_t net = leaf_growth_rate + root_growth_rate;
    /* Senescence returns tissue; it does not refund construction carbon. */
    if (net <= REAL(0.0))
        return REAL(0.0);
    return growth_stoichiometry * net;
}

real_t transport_cost(
    real_t sucrose_consumption_transport,
    real_t respiration_frequency,
    real_t sucrose,
    real_t sucrose_loading_frequency,
    real_t night_efficiency_starch)
{
    return sucrose_consumption_transport * ((respiration_frequency * sucrose) + (sucrose_loading_frequency * night_efficiency_starch));
}

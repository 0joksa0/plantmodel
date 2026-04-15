#include "model/model.h"

real_t limitation_of_photosynthetic_rate(
    real_t starch,
    real_t feedback_on_photosynthesis,
    real_t max_starch)
{
    return feedback_on_photosynthesis + ((REAL(1.0) - feedback_on_photosynthesis) * ((max_starch - starch) / max_starch));
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
        input->core.light,
        starch_partition_coeff,
        input->core.min_sucrose,
        input->core.max_sucrose,
        input->core.sucrose,
        input->core.lambda_sdr,
        input->core.lambda_sdi,
        input->core.lambda_sni);
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
    return starch_production(input->core.photosynthesis, input->core.starch_partition_coeff, input->core.starch_degradation_rate);
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
    real_t night_efficiency_starch)
{
    return ((REAL(1.0) - sucrose_partition_coeff) * photosynthesis) + starch_degradation_rate - uptake_cost - transport_cost - (respiration_frequency * sucrose) - (sucrose_loading_frequency * night_efficiency_starch);
}

real_t sucrose_production_f(
    real_t t,
    real_t sucrose,
    void* params)
{
    Input* in = (Input*)params;
    (void)t;

    return sucrose_production(in->core.starch_partition_coeff,
        in->core.photosynthesis,
        in->core.starch_degradation_rate,
        in->core.uptake_cost,
        in->core.transport_cost,
        in->core.respiration_frequency,
        sucrose,
        in->core.sucrose_loading_frequency,
        in->core.night_efficiency_starch);
}

real_t night_efficiency_starch(
    real_t sucrose,
    real_t max_sucrose,
    real_t lambda_g,
    real_t light,
    real_t starch_partition_coeff)
{

    return RMIN(sucrose, max_sucrose) * (lambda_g + ((REAL(1.0) - lambda_g) * (REAL(1.0) - light + (light * (REAL(1.0) - starch_partition_coeff)))));
}

real_t uptake_cost(
    real_t nitrogen_uptake_sucrose_consumption,
    real_t nitrogen_uptake,
    real_t phosphorus_uptake_sucrose_consumption,
    real_t phosphorus_uptake)
{
    return (nitrogen_uptake_sucrose_consumption * nitrogen_uptake) + (phosphorus_uptake_sucrose_consumption * phosphorus_uptake);
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

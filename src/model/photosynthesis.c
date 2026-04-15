#include "model/model.h"

#include <math.h>
#include <stdio.h>

static FILE *g_photo_fp = NULL;

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

real_t farquhar_photosynthesis(Input *input)
{
    if (input->photo.light <= REAL(0.0) || input->growth.leaf_biomass < input->photo.min_leaf_biomass)
    {
        return REAL(0.0);
    }
    if (input->photo.light_PAR < REAL(20.0))
    {
        return REAL(0.0);
    }

    real_t VPD = calculate_vapor_pressure_deficit(input->gas_exchange.leaf_temperature, input->gas_exchange.relative_humidity);
    input->gas_exchange.stomatal_conductance = calculate_stomatal_conductance_jarvis(
        input->photo.light_PAR,
        VPD,
        input->gas_exchange.leaf_temperature,
        input->gas_exchange.ambient_CO2_concentration,
        input->gas_exchange.leaf_water_potential);

    real_t An = iterate_ci(input, 1000, REAL(0.0001));

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

    return RMAX(REAL(0.0), adjusted_PH);
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
    const real_t base_scale = REAL(1.10);
    const real_t ref_specific_leaf_area = REAL(0.025);
    real_t safe_sla = RMAX(specific_leaf_area, REAL(1e-6));
    real_t sla_scale = safe_sla / ref_specific_leaf_area;

    (void)leaf_biomass;
    return net_photosynthesis_rate * base_scale * sla_scale;
}

real_t intercellular_CO2(
    real_t ambient_CO2_concentration,
    real_t net_photosynthesis_rate,
    real_t stomatal_conductance)
{
    return ambient_CO2_concentration - (net_photosynthesis_rate / stomatal_conductance);
}

real_t iterate_ci(Input *input, int max_iter, real_t epsilon)
{
    real_t Ci = input->gas_exchange.ambient_CO2_concentration * REAL(0.7);
    real_t prev_An = 0.0;
    real_t An = 0.0;
    real_t gs = RMAX(input->gas_exchange.stomatal_conductance, REAL(1e-6));

    for (int i = 0; i < max_iter; i++)
    {
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

        real_t new_Ci = input->gas_exchange.ambient_CO2_concentration - An / gs;
        new_Ci = RMAX(REAL(0.0), new_Ci);

        if (RABS(new_Ci - Ci) < epsilon && RABS(An - prev_An) < epsilon)
        {
            break;
        }

        Ci = new_Ci;
        prev_An = An;
    }

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

    const real_t T0 = REAL(25.0);
    const real_t T1 = REAL(5.0);
    const real_t Th = REAL(45.0);
    const real_t b4 = (Th - T0) / (Th - T1);
    const real_t b3 = REAL(1.0) / ((T0 - T1) * RPOW((Th - T1), b4));
    real_t f_T = b3 * (leaf_temperature_celsius - T1) * RPOW((Th - leaf_temperature_celsius), b4);
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

    return f_PAR * f_VPD * f_T * f_CO2 * f_water;
}

real_t calculate_vapor_pressure_deficit(real_t temperature_celsius, real_t relative_humidity_percent)
{
    real_t saturation_vapor_pressure = REAL(0.61078) * REXP((REAL(17.27) * temperature_celsius) / (temperature_celsius + REAL(237.3)));
    real_t actual_vapor_pressure = saturation_vapor_pressure * (relative_humidity_percent / REAL(100.0));
    return saturation_vapor_pressure - actual_vapor_pressure;
}

void update_light_conditions(Input *input, real_t current_hour)
{
    const real_t max_PAR = 1200.0;

    real_t sunrise = 0.0;
    real_t sunset = sunrise + input->photo.photoperiod;

    if (current_hour >= sunrise && current_hour < sunset)
    {
        input->photo.light = 1.0;
        real_t day_fraction = (current_hour - sunrise) / (sunset - sunrise);
        input->photo.light_PAR = max_PAR * RSIN(day_fraction * M_PI);
    }
    else
    {
        input->photo.light = 0.0;
        input->photo.light_PAR = 0.0;
    }
}

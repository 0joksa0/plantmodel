#include "export/export.h"
#include "gui/plot.h"
#include "model/input.h"
#include "model/model.h"
#include "simulation/runge_cutta.h"
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#define SIMULATION_DT 0.1f

int total_steps = 0;
float* sucrose = NULL;
float* starch = NULL;
float* ph = NULL;
float* partition = NULL;
float eps = 1e-6f;
bool trigger = false;
float rgr(float fw_now, float fw_start, float hours_passed)
{
    return (logf(fw_now) - logf(fw_start)) / hours_passed;
}

/* poziv na kraju simulacije */

void simulate_step(float current_time, int step, Input* input)
{
    // --- Osvetljenje na osnovu trenutnog vremena ---
    float hour_of_day = fmodf(current_time, 24.0f);
    input->light = (hour_of_day >= 6.0f && hour_of_day < 18.0f) ? 1 : 0;
    if (input->light == 0) {
        if (!trigger) {
            input->starch_night_start = input->starch;
        }
        trigger = true;
    } else
        trigger = false;

    // --- Ograničenje fotosinteze ---
    input->limitation_of_photosyntetic_rate = limitation_of_photosyntethic_rate(
        input->starch, input->feedback_on_photosynthesis, input->max_starch);

    input->nitrogen_saturation = nitrogen_saturation(
        input->nitrogen, input->min_nitrogen_photosynthesis);

    input->phosphorus_saturation = phosphorus_saturation(
        input->min_nitrogen_photosynthesis, input->optimal_stechiometric_ratio, input->phosphorus);

    // --- Fotosinteza ---
    input->photosynthesis = photosynthesis(
        input->light,
        input->limitation_of_photosyntetic_rate,
        input->max_photosyntetic_rate,
        input->nitrogen_saturation,
        input->phosphorus_saturation,
        input->leaf_biomass,
        input->min_leaf_biomass);

    // --- Degradacija skroba ---
    input->starch_degradation_rate = starch_degradation(
        input->max_starch_degradation_rate,
        input->max_sucrose,
        input->starch_night_start,
        input->min_starch,
        input->sucrose,
        input->light,
        hour_of_day,
        input->photoperiod);

    // --- Podela saharoze/skroba ---
    input->starch_partition_coeff = runge_kutta_4(input->starch_partition_coeff, current_time, SIMULATION_DT, starch_sucrose_partition_f, input);

    // --- Update biomasa i resursa uz RK4 ---
    // 7. Proizvodnja skroba i saharoze
    input->starch = runge_kutta_4(input->starch, current_time, SIMULATION_DT, starch_production_f, input);
    input->sucrose = runge_kutta_4(input->sucrose, current_time, SIMULATION_DT, sucrose_production_f, input);

    // 8-10. Uptake azota i fosfora
    input->nitrogen_nutrient_uptake = nitrogen_nutrient_uptake(input->max_nitrogen_uptake, input->nitrogen_soil_content, input->Michaelis_Menten_constant_nitrogen);
    input->phosphorus_nutrient_uptake = phosphorus_nutrient_uptake(input->max_phosphorus_uptake, input->phosphorus_soil_content, input->Michaelis_Menten_constant_phosphorus);

    input->pot_nitrogen_uptake = pot_nitrogen_uptake(input->nitrogen_nutrient_uptake, input->nitrogen_affinity, input->root_biomass, input->total_biomass);
    input->pot_phosphorus_uptake = pot_phosphorus_uptake(input->phosphorus_nutrient_uptake, input->phosphorus_affinity, input->root_biomass, input->total_biomass);

    input->nitrogen_uptake = nitrogen_uptake(input->pot_nitrogen_uptake, input->sucrose, input->nitrogen_uptake_sucrose_consumption);
    input->phosphorus_uptake = phosphorus_uptake(input->pot_phosphorus_uptake, input->sucrose, input->phosphorus_uptake_sucrose_consumption);

    // 11. Korekcije azota i fosfora
    input->nitrogen = runge_kutta_4(input->nitrogen, current_time, SIMULATION_DT, nitrogen_content_f, input);
    input->phosphorus = runge_kutta_4(input->phosphorus, current_time, SIMULATION_DT, phosphorus_content_f, input);

    // 12. Troškovi uptake-a
    input->uptake_cost = uptake_cost(input->nitrogen_uptake_sucrose_consumption, input->nitrogen_uptake, input->phosphorus_uptake_sucrose_consumption, input->phosphorus_uptake);

    // 13. Troškovi transporta
    input->transport_cost = transport_cost(input->sucrose_consumption_transport, input->respiration_frequency, input->sucrose, input->sucrose_loading_frequency, input->night_efficiency_starch);

    input->nitrogen_cost = nitrogen_cost(input->respiration_frequency, input->sucrose, input->starch, input->sucrose_loading_frequency, input->night_efficiency_starch, input->assimilation_cost_nitrogen, input->leaf_biomass, input->total_biomass, input->photosynthesis, input->max_photosyntetic_rate, input->min_nitrogen_photosynthesis);
    input->phosphorus_cost = phosphorus_cost(input->respiration_frequency, input->sucrose, input->starch, input->sucrose_loading_frequency, input->night_efficiency_starch, input->assimilation_cost_phorphorus, input->leaf_biomass, input->total_biomass, input->photosynthesis, input->max_photosyntetic_rate, input->min_phosphorus_photosynthesis);

    // 15. Korekcije afiniteta
    input->nitrogen_affinity = runge_kutta_4(input->nitrogen_affinity, current_time, SIMULATION_DT, nitrogen_affinity_f, input);
    input->phosphorus_affinity = runge_kutta_4(input->phosphorus_affinity, current_time, SIMULATION_DT, phosphorus_affinity_f, input);

    // Clamp za afinitet

    // 16. Night efficiency
    input->night_efficiency_starch = night_efficieny_starch(input->sucrose, input->max_sucrose, input->lambda_g, input->light, input->starch_partition_coeff);

    // 17. Rast lišća
    input->leaf_biomass = runge_kutta_4(input->leaf_biomass, current_time, SIMULATION_DT, leaf_growth_f, input);

    // 18. Rast korena
    input->root_biomass = runge_kutta_4(input->root_biomass, current_time, SIMULATION_DT, root_growth_f, input);

    // 19. Ukupna biomasa
    input->total_biomass = input->leaf_biomass + input->root_biomass;
    input->max_starch = max_starch(input->max_starch_degradation_rate, input->photoperiod);
    input->stochiometric_signal = stochiometric_signal(input->optimal_stechiometric_ratio, input->nitrogen, input->phosphorus);

    // 20. Korekcija alokacije sukroze
    input->sucrose_root_allocation = runge_kutta_4(input->sucrose_root_allocation, current_time, SIMULATION_DT, sucrose_root_allocation_f, input);

    // Clamp za sukrozu root alokaciju

    // 21. Korekcija stoichiometric signal

    // 17. Loguj stanje
    log_simulation_step(current_time, input);
}

void simulate_days(int days)
{
    total_steps = (int)(days * 24.0f / SIMULATION_DT);
    float current_time = 0.0f;
    sucrose = malloc(total_steps * sizeof(float));
    starch = malloc(total_steps * sizeof(float));

    ph = malloc(total_steps * sizeof(float));

    partition = malloc(total_steps * sizeof(float));

    Input input = generate_input();

    float start_leafe_biomass = input.leaf_biomass;

    for (int step = 0; step < total_steps; step++) {
        simulate_step(current_time, step, &input);
        current_time += SIMULATION_DT;
        sucrose[step] = input.sucrose;
        starch[step] = input.starch;
        ph[step] = input.photosynthesis;
        partition[step] = input.total_biomass;
        // printf("%d\n", step);
    }

    main_thread();
    /* simulation.c, posle simulacije */
    printf("t_end - t0  = %.0f h\n", days * 24.0f); //  EXPECT 696
    printf("currentTime = %d h\n", total_steps); //  Šta zapravo dobiješ?

    float RGR_total = rgr(input.leaf_biomass, start_leafe_biomass, current_time);
    printf("leaf FW in t0: %f, leaf FW u t: %f, t0:%f, t:%f, dt:%f", start_leafe_biomass, input.leaf_biomass, 0.0f, current_time, SIMULATION_DT);
    printf("\nRGR: %f  RGR_F: %f \n", RGR_total / 0.018, RGR_total * 24.0f);
}

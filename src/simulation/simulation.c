#include "export/export.h"
#include "gui/plot.h"
#include "model/input.h"
#include "model/model.h"
#include <solver.h>

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#define SIMULATION_DT REAL(1.0)
#define PHOTO_PERIOD REAL(8.0)

SolverType solver_type = SOLVER_ODE45;

real_t solver_tolerance = REAL(1e-8);

int total_steps = 0;
real_t* sucrose = NULL;
real_t* starch = NULL;
real_t* ph = NULL;
real_t* partition = NULL;
real_t eps = REAL_EPSILON;

bool trigger = false;

real_t rgr(real_t fw_now, real_t fw_start, real_t hours_passed)
{
    return (log(fw_now / fw_start)) / hours_passed;
}

void simulate_step(real_t current_time, int step, Input* input)
{
    real_t hour_of_day = fmod(current_time, REAL(24.0));
    input->light = (hour_of_day >= REAL(0.0) && hour_of_day <= PHOTO_PERIOD) ? 1 : 0;
    if (input->light == 0) {
        if (!trigger) {
            input->starch_night_start = input->starch;
        }
        trigger = true;

    } else
        trigger = false;

    input->min_nitrogen = min_nitrogen(input->respiration_frequency, input->max_sucrose, input->max_starch, input->sucrose_loading_frequency, input->assimilation_cost_nitrogen, input->min_nitrogen_photosynthesis, input->nutrient_conversion_parameter, input->photoperiod);
    input->max_nitrogen = max_nitrogen(input->min_nitrogen);

    input->min_phosphorus = min_phosphorus(
        input->respiration_frequency, input->max_sucrose, input->max_starch, input->sucrose_loading_frequency, input->assimilation_cost_phorphorus, input->min_phosphorus_photosynthesis, input->nutrient_conversion_parameter, input->photoperiod);
    input->max_phosphorus = max_phosphorus(input->max_nitrogen, input->optimal_stechiometric_ratio);

    input->nitrogen_saturation = nitrogen_saturation(
        input->nitrogen, input->min_nitrogen_photosynthesis);
    input->phosphorus_saturation = phosphorus_saturation(
        input->min_nitrogen_photosynthesis, input->optimal_stechiometric_ratio, input->phosphorus);

    input->limitation_of_photosyntetic_rate = limitation_of_photosyntethic_rate(
        input->starch, input->feedback_on_photosynthesis, input->max_starch);

    // input->photosynthesis = photosynthesis(
    //     input->light,
    //     input->limitation_of_photosyntetic_rate,
    //     input->max_photosyntetic_rate,
    //     input->nitrogen_saturation,
    //     input->phosphorus_saturation,
    //     input->leaf_biomass,
    //     input->min_leaf_biomass);
    input->photosynthesis = farquhar_photosynthesis(input);

    input->night_efficiency_starch = night_efficieny_starch(input->sucrose, input->max_sucrose, input->lambda_g, input->light, input->starch_partition_coeff);

    input->uptake_cost = uptake_cost(input->nitrogen_uptake_sucrose_consumption, input->nitrogen_uptake, input->phosphorus_uptake_sucrose_consumption, input->phosphorus_uptake);
    input->transport_cost = transport_cost(input->sucrose_consumption_transport, input->respiration_frequency, input->sucrose, input->sucrose_loading_frequency, input->night_efficiency_starch);
    input->nitrogen_cost = nitrogen_cost(input->respiration_frequency, input->sucrose, input->starch, input->sucrose_loading_frequency, input->night_efficiency_starch, input->assimilation_cost_nitrogen, input->leaf_biomass, input->total_biomass, input->photosynthesis, input->max_photosyntetic_rate, input->min_nitrogen_photosynthesis);
    input->phosphorus_cost = phosphorus_cost(input->respiration_frequency, input->sucrose, input->starch, input->sucrose_loading_frequency, input->night_efficiency_starch, input->assimilation_cost_phorphorus, input->leaf_biomass, input->total_biomass, input->photosynthesis, input->max_photosyntetic_rate, input->min_phosphorus_photosynthesis);

    input->starch_degradation_rate = starch_degradation(
        input->max_starch_degradation_rate,
        input->max_sucrose,
        input->starch_night_start,
        input->min_starch,
        input->sucrose,
        input->light,
        hour_of_day,
        input->photoperiod);



    input->starch_partition_coeff = solve_step(solver_type, solver_tolerance, input->starch_partition_coeff, current_time, SIMULATION_DT, starch_sucrose_partition_f, input);

    input->starch = solve_step(solver_type, solver_tolerance, input->starch, current_time, SIMULATION_DT, starch_production_f, input);
    input->sucrose = solve_step(solver_type, solver_tolerance, input->sucrose, current_time, SIMULATION_DT, sucrose_production_f, input);

    input->nitrogen_affinity = solve_step(solver_type, solver_tolerance, input->nitrogen_affinity, current_time, SIMULATION_DT, nitrogen_affinity_f, input);
    input->phosphorus_affinity = solve_step(solver_type, solver_tolerance, input->phosphorus_affinity, current_time, SIMULATION_DT, phosphorus_affinity_f, input);

    input->nitrogen_nutrient_uptake = nitrogen_nutrient_uptake(input->max_nitrogen_uptake, input->nitrogen_soil_content, input->Michaelis_Menten_constant_nitrogen);
    input->phosphorus_nutrient_uptake = phosphorus_nutrient_uptake(input->max_phosphorus_uptake, input->phosphorus_soil_content, input->Michaelis_Menten_constant_phosphorus);

    input->pot_nitrogen_uptake = pot_nitrogen_uptake(input->nitrogen_nutrient_uptake, input->nitrogen_affinity, input->root_biomass, input->total_biomass);
    input->pot_phosphorus_uptake = pot_phosphorus_uptake(input->phosphorus_nutrient_uptake, input->phosphorus_affinity, input->root_biomass, input->total_biomass);

    input->nitrogen_uptake = nitrogen_uptake(input->pot_nitrogen_uptake, input->sucrose, input->nitrogen_uptake_sucrose_consumption);
    input->phosphorus_uptake = phosphorus_uptake(input->pot_phosphorus_uptake, input->sucrose, input->phosphorus_uptake_sucrose_consumption);

    input->nitrogen = solve_step(solver_type, solver_tolerance, input->nitrogen, current_time, SIMULATION_DT, nitrogen_content_f, input);
    input->phosphorus = solve_step(solver_type, solver_tolerance, input->phosphorus, current_time, SIMULATION_DT, phosphorus_content_f, input);

    // printf("n_aff:%f,n_nut_uptake: %f, n_pot_uptake: %f, n_uptake: %f, n: %f\n", input->nitrogen_affinity, input->nitrogen_nutrient_uptake, input->pot_nitrogen_uptake, input->nitrogen_uptake, input->nitrogen);
    input->max_starch = max_starch(input->max_starch_degradation_rate, input->photoperiod);
    input->stochiometric_signal = stochiometric_signal(input->optimal_stechiometric_ratio, input->nitrogen, input->phosphorus);

    input->sucrose_root_allocation = solve_step(solver_type, solver_tolerance, input->sucrose_root_allocation, current_time, SIMULATION_DT, sucrose_root_allocation_f, input);

    input->leaf_biomass = solve_step(solver_type, solver_tolerance, input->leaf_biomass, current_time, SIMULATION_DT, leaf_growth_f, input);

    input->root_biomass = solve_step(solver_type, solver_tolerance, input->root_biomass, current_time, SIMULATION_DT, root_growth_f, input);

    input->total_biomass = input->leaf_biomass + input->root_biomass;
    log_simulation_step(current_time, input);
}

void simulate_days(int days, Input* input)
{
    total_steps = (int)(days * REAL(24.0) / SIMULATION_DT);
    real_t current_time = REAL(0.0);
    sucrose = malloc(total_steps * sizeof(real_t));
    starch = malloc(total_steps * sizeof(real_t));

    ph = malloc(total_steps * sizeof(real_t));

    partition = malloc(total_steps * sizeof(real_t));

    real_t start_biomass = input->leaf_biomass;
    input->photoperiod = PHOTO_PERIOD;

    for (int step = 0; step < total_steps; step++) {
        simulate_step(current_time, step, input);
        current_time += SIMULATION_DT;
        sucrose[step] = input->sucrose;
        starch[step] = input->starch;
        ph[step] = input->photosynthesis;
        partition[step] = input->total_biomass;
        printf("\n%d\n", step);
    }

    main_thread();
    /* simulation.c, posle simulacije */
    printf("t_end - t0  = %.0f h\n", days * REAL(24.0)); //  EXPECT 696
    printf("currentTime = %d h\n", total_steps); //  Šta zapravo dobiješ?

    real_t RGR_total = rgr(input->leaf_biomass, start_biomass, days);
    printf("leaf FW in t0: %f, leaf FW u t: %f, t0:%f, t:%f, dt:%f", start_biomass, input->leaf_biomass, REAL(0.0), current_time, SIMULATION_DT);
    printf("\nRGR: %f  RGR_F: %f \n", RGR_total, RGR_total * REAL(24.0));

    printf("total_biomass: %f", input->total_biomass);
}

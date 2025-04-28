#include "../include/model.h"
#include "../include/runge_cutta.h"
#include <math.h>
#include <raylib.h>
#include <stdio.h>

#define SIMULATION_DT 0.1f // novi vremenski korak 0.1h (6 minuta)

float eps = 1e-6f;

void DrawAxisGrid(int scrollX, int scrollY, float zoom)
{
    int screenWidth = GetScreenWidth();
    int screenHeight = GetScreenHeight();

    // X-axis
    DrawLine(0, 500 - scrollY, screenWidth, 500 - scrollY, LIGHTGRAY);
    for (int i = -scrollX; i < screenWidth / zoom + scrollX; i += 24) {
        int worldX = i;
        int screenX = 20 + worldX * zoom - scrollX;
        DrawLine(screenX, 495 - scrollY, screenX, 505 - scrollY, GRAY);
        DrawText(TextFormat("%d", worldX), screenX - 10, 510 - scrollY, 10, DARKGRAY);
    }

    // Y-axis
    DrawLine(20 - scrollX, 0, 20 - scrollX, screenHeight, LIGHTGRAY);
    for (int j = -scrollY; j < screenHeight / (20 * zoom) + scrollY; j += 5) {
        int worldY = j;
        int screenY = 500 - worldY * 20 * zoom - scrollY;
        DrawLine(15 - scrollX, screenY, 25 - scrollX, screenY, GRAY);
        DrawText(TextFormat("%d", worldY), 0, screenY - 5, 10, DARKGRAY);
    }
}
void LogSimulationStep(int t, Input* input)
{
    t--;
    FILE* log = fopen("sim_log.csv", t == 0 ? "w" : "a");
    if (log) {
        if (t == 0) {
            fprintf(log,
                "Hour,Light,Gamma,Photosynthesis,"
                "Sucrose,Starch,StarchDegr,StarchProd,SucroseProd,"
                "Nitrogen,Phosphorus,N_Aff,P_Aff,"
                "RespCost,UptakeCost,TransportCost,"
                "LeafBiomass,RootBiomass,TotalBiomass,Photoperiod\n");
        }

        fprintf(log,
            "%d,%.0f,%.4f,%.4f,"
            "%.4f,%.4f,%.4f,%.4f,%.4f,"
            "%.4f,%.4f,%.4f,%.4f,"
            "%.4f,%.4f,%.4f,"
            "%.4f,%.4f,%.4f,%.2f\n",
            t, input->light,
            input->starch_partition_coeff, input->photosynthesis,
            input->sucrose, input->starch,
            input->starch_degradation_rate, // ovo je degradacija
            input->starch_partition_coeff * input->photosynthesis - input->starch_degradation_rate, // starch_prod
            (1.0f - input->starch_partition_coeff) * input->photosynthesis + input->starch_degradation_rate - input->uptake_cost - input->transport_cost - (input->respiration_frequency * input->sucrose) - (input->sucrose_loading_frequency * input->night_efficiency_starch), // sucrose_prod
            input->nitrogen, input->phosphorus,
            input->nitrogen_affinity, input->phosphorus_affinity,
            input->night_efficiency_starch, input->uptake_cost, input->transport_cost,
            input->leaf_biomass, input->root_biomass, input->total_biomass,
            input->photoperiod);

        fclose(log);
    }
}

void simulate_step(float current_time, int step, Input* input)
{
    // --- Osvetljenje na osnovu trenutnog vremena ---
    float hour_of_day = fmodf(current_time, 24.0f);
    printf("\ncurrent time: %f, hour of day: %f \n", current_time, hour_of_day);
    input->light = (hour_of_day >= 6.0f && hour_of_day < 18.0f) ? 1 : 0;
    if (input->light == 0)
        input->starch_night_start = input->starch;
    printf("light : %f", input->light);

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
    input->stochiometric_signal = stochiometric_signal(input->optimal_stechiometric_ratio, input->nitrogen, input->phosphorus);

    // 20. Korekcija alokacije sukroze
    input->sucrose_root_allocation = runge_kutta_4(input->sucrose_root_allocation, current_time, SIMULATION_DT, sucrose_root_allocation_f, input);

    // Clamp za sukrozu root alokaciju

    // 21. Korekcija stoichiometric signal

    // 17. Loguj stanje
    LogSimulationStep(step, input);
}

void simulate_days(int days)
{
    int total_steps = (int)(days * 24.0f / SIMULATION_DT);
    float current_time = 0.0f;
    float sucroseValues[10000] = { 0 };
    float starchValues[10000] = { 0 };
    float phValues[10000] = { 0 };
    float partition[10000] = { 0 };

    Input input = {
        .light = 0.0f,
        .limitation_of_photosyntetic_rate = 0.5f,
        .max_photosyntetic_rate = 12.7f,
        .nitrogen_saturation = 1.0f,
        .phosphorus_saturation = 1.0f,
        .leaf_biomass = 0.6f,
        .root_biomass = 1.0f,
        .total_biomass = 2.0f,
        .min_leaf_biomass = 0.5f,
        .feedback_on_photosynthesis = 0.5f,
        .max_starch = 72.0f,
        .max_starch_degradation_rate = 6.0f,
        .photoperiod = 12.0f,
        .optimal_stechiometric_ratio = 10.0f,
        .starch_partition_coeff = 0.5f,
        .starch_degradation_rate = 0.0f,
        .uptake_cost = 0.0f,
        .transport_cost = 0.0f,
        .respiration_frequency = 0.79f,
        .sucrose_loading_frequency = 1.98f,
        .night_efficiency_starch = 0.5f,
        .max_sucrose = 2.0f,
        .min_sucrose = 1.3f,
        .starch_night_start = 2.0f,
        .min_starch = 0.15f,
        .lambda_sdr = 0.25f,
        .lambda_sdi = 0.10f,
        .lambda_sni = 0.08f,
        .lambda_g = 0.1f,
        .lambda_sb = 0.00587f,
        .sucrose = 1.5f,
        .starch = 2.0f,
        .nitrogen = 6.5f,
        .phosphorus = 0.65f,
        .nitrogen_uptake = 0.0f,
        .phosphorus_uptake = 0.0f,
        .assimilation_cost_nitrogen = 0.04758f,
        .lambda_csn = 0.0267f,
        .min_nitrogen_photosynthesis = 0.75f,
        .min_phosphorus_photosynthesis = 0.075f,
        .nutrient_conversion_parameter = 0.066f,
        .assimilation_cost_phorphorus = 0.004758f,
        .pot_nitrogen_uptake = 0.0f,
        .pot_phosphorus_uptake = 0.0f,
        .nitrogen_nutrient_uptake = 0.0f,
        .phosphorus_nutrient_uptake = 0.0f,
        .nitrogen_affinity = 0.5f,
        .phosphorus_affinity = 0.5f,
        .nitrogen_soil_content = 10.0f,
        .phosphorus_soil_content = 0.15f,
        .nitrogen_uptake_sucrose_consumption = 0.11f,
        .phosphorus_uptake_sucrose_consumption = 0.11f,
        .max_nitrogen_uptake = 6.44f,
        .max_phosphorus_uptake = 0.4f,
        .Michaelis_Menten_constant_nitrogen = 0.125f,
        .Michaelis_Menten_constant_phosphorus = 0.006736f,
        .lambda_k = 1.0f,
        .max_nitrogen = 23.0f,
        .min_nitrogen = 5.75f,
        .max_phosphorus = 2.3f,
        .min_phosphorus = 0.575f,
        .nitrogen_cost = 0.0f,
        .phosphorus_cost = 0.0f,
        .sucrose_consumption_transport = 0.0035f,
        .sucrose_root_allocation = 0.5f,
        .stochiometric_signal = 0.5f,
        .leaf_deathrate = 0.000005f,
        .root_deathrate = 0.000005f,
        .leaf_competitive_rate = 0.000005f,
        .root_competitive_rate = 0.000005f
    };

    InitWindow(GetMonitorWidth(GetCurrentMonitor()), GetMonitorHeight(GetCurrentMonitor()), "Plant Simulation");
    SetTargetFPS(60);
    int scrollX = 0;
    int scrollY = 0;
    int zoom = 5.0f;
    for (int step = 0; step < total_steps; step++) {
        simulate_step(current_time, step, &input);
        current_time += SIMULATION_DT;
        sucroseValues[step] = input.sucrose;
        starchValues[step] = input.starch;
        phValues[step] = input.photosynthesis;
        partition[step] = input.total_biomass;
        printf("%d\n", step);
    }

    while (!WindowShouldClose()) {
        if (IsKeyDown(KEY_RIGHT))
            scrollX += 5;
        if (IsKeyDown(KEY_LEFT))
            scrollX -= 5;
        if (IsKeyDown(KEY_DOWN))
            scrollY += 5;
        if (IsKeyDown(KEY_UP))
            scrollY -= 5;
        if (IsKeyPressed(KEY_A))
            zoom *= 2.0f;
        if (IsKeyPressed(KEY_B))
            zoom /= 1.1f;

        BeginDrawing();
        ClearBackground(RAYWHITE);
        DrawAxisGrid(scrollX, scrollY, zoom);
        for (int i = 1; i < total_steps; ++i) {
            DrawLine(2 + (i - 1 - scrollX) * zoom, 50 - sucroseValues[i - 1] * 20 * zoom - scrollY, 2 + (i - scrollX) * zoom, 50 - sucroseValues[i] * 20 * zoom - scrollY, RED);
            DrawLine(2 + (i - 1 - scrollX) * zoom, 50 - starchValues[i - 1] * 20 * zoom - scrollY, 2 + (i - scrollX) * zoom, 50 - starchValues[i] * 20 * zoom - scrollY, BLUE);
            DrawLine(2 + (i - 1 - scrollX) * zoom, 50 - phValues[i - 1] * 20 * zoom - scrollY, 2 + (i - scrollX) * zoom, 50 - phValues[i] * 20 * zoom - scrollY, GREEN);
            DrawLine(2 + (i - 1 - scrollX) * zoom, 50 - partition[i - 1] * 20 * zoom - scrollY, 2 + (i - scrollX) * zoom, 50 - partition[i] * 20 * zoom - scrollY, BLACK);
        }

        DrawText("Sucrose", 20, 20, 20, RED);
        DrawText("Starch", 20, 40, 20, BLUE);
        DrawText("Photosynthesis", 20, 60, 20, GREEN);
        DrawText("Biomass", 20, 80, 20, BLACK);
        DrawText(TextFormat("Zoom: %.2fx", zoom), 1000, 20, 20, DARKGRAY);

        EndDrawing();
    }
    CloseWindow();
}

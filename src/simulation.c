#include "../include/model.h"
// simulation.c
#include <math.h>
#include <raylib.h>
#include <stdio.h>

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

void simulate_step(int t, Input* input)
{
    float dt = 1.0f; // Simulacioni korak: 1 sat
    input->light = (t % 24 >= 6 && t % 24 < 18) ? 1.0f : 0.0f;

    // 2. Ograničenje biomase (nema negativne biomase)
    // 3. Limitacija fotosinteze
    input->limitation_of_photosyntetic_rate = limitation_of_photosyntethic_rate(
        input->starch, input->feedback_on_photosynthesis, input->max_starch);
    
    printf("lim_phot_rate: %f", input->limitation_of_photosyntetic_rate);

    input->nitrogen_saturation = nitrogen_saturation(
        input->nitrogen, input->min_nitrogen_photosynthesis);

    printf("nitrogen_sat: %f", input->nitrogen_saturation);
    input->phosphorus_saturation = phosphorus_saturation(
        input->min_nitrogen_photosynthesis, input->optimal_stechiometric_ratio, input->phosphorus);

    printf("phosph_sat: %f", input->phosphorus_saturation);
    // 4. Fotosinteza
    float ph = photosynthesis(
        input->light,
        input->limitation_of_photosyntetic_rate,
        input->max_photosyntetic_rate,
        input->nitrogen_saturation,
        input->phosphorus_saturation,
        input->leaf_biomass,
        input->min_leaf_biomass);
    printf("fotosinteza: %f", ph);
    input->photosynthesis = ph;
    // 5. Degradacija skroba
    input->starch_degradation_rate = starch_degradation(
        input->max_starch_degradation_rate,
        input->max_sucrose,
        input->starch_night_start,
        input->min_starch,
        input->sucrose,
        input->light,
        (float)(t % 24),
        input->photoperiod);

    printf("starchdeg: %f", input->starch_degradation_rate);
    // 6. Podela skroba/saharoze (particija)
    input->starch_partition_coeff = starch_sucrose_partition(
        input->light,
        input->starch_partition_coeff,
        input->min_sucrose,
        input->max_sucrose,
        input->sucrose,
        input->lambda_sdr,
        input->lambda_sdi,
        input->lambda_sni);

    printf("starch suc part: %f", input->starch_partition_coeff);
    // 7. Proizvodnja skroba i saharoze
    float starch_prod = starch_production(input->photosynthesis, input->starch_partition_coeff, input->starch_degradation_rate);
    printf("starch prod: %f", starch_prod);
    float sucrose_prod = sucrose_production(input->starch_partition_coeff, input->photosynthesis, input->starch_degradation_rate, input->uptake_cost, input->transport_cost, input->respiration_frequency, input->sucrose, input->sucrose_loading_frequency, input->night_efficiency_starch);
    printf("sucrose: %f", sucrose_prod);

    // 8. Uptake i troškovi
    input->nitrogen_nutrient_uptake = nitrogen_nutrient_uptake(
        input->max_nitrogen_uptake,
        input->nitrogen_soil_content,
        input->Michaelis_Menten_constant_nitrogen);
    input->phosphorus_nutrient_uptake = phosphorus_nutrient_uptake(
        input->max_phosphorus_uptake,
        input->phosphorus_soil_content,
        input->Michaelis_Menten_constant_phosphorus);

    input->pot_nitrogen_uptake = pot_nitrogen_uptake(
        input->nitrogen_nutrient_uptake,
        input->nitrogen_affinity,
        input->root_biomass,
        input->total_biomass);

    input->pot_phosphorus_uptake = pot_phosphorus_uptake(
        input->phosphorus_nutrient_uptake,
        input->phosphorus_affinity,
        input->root_biomass,
        input->total_biomass);

    // 9. Korekcije azota i fosfora
    input->nitrogen_uptake = nitrogen_uptake(
        input->pot_nitrogen_uptake,
        input->sucrose,
        input->nitrogen_uptake_sucrose_consumption);
    input->phosphorus_uptake = phosphorus_uptake(
        input->pot_nitrogen_uptake,
        input->sucrose,
        input->phosphorus_uptake_sucrose_consumption);

    // 10. Ažuriranje biomasa i resursa
    input->sucrose = input->sucrose + sucrose_prod;
    input->starch = input->starch + starch_prod;
    input->nitrogen = input->nitrogen + input->nitrogen_uptake;
    input->phosphorus = input->phosphorus + input->phosphorus_uptake;

    // 11. Troškovi transporta i disanja
    input->uptake_cost = uptake_cost(
        input->nitrogen_uptake_sucrose_consumption,
        input->nitrogen_uptake,
        input->phosphorus_uptake_sucrose_consumption,
        input->phosphorus_uptake);
    input->nitrogen_cost = nitrogen_cost(
        input->respiration_frequency,
        input->sucrose,
        input->starch,
        input->sucrose_loading_frequency,
        input->night_efficiency_starch,
        input->assimilation_cost_nitrogen,
        input->leaf_biomass,
        input->total_biomass,
        input->photosynthesis,
        input->max_photosyntetic_rate,
        input->min_nitrogen_photosynthesis);

    input->phosphorus_cost = phosphorus_cost(
        input->respiration_frequency,
        input->sucrose,
        input->starch,
        input->sucrose_loading_frequency,
        input->night_efficiency_starch,
        input->assimilation_cost_phorphorus,
        input->leaf_biomass,
        input->total_biomass,
        input->photosynthesis,
        input->max_photosyntetic_rate,
        input->min_phosphorus_photosynthesis);

    // input->nitrogen_affinity = nitrogen_affinity(
    //  input->nitrogen_affinity,
    //  input->nitrogen_uptake,
    //  input->nitrogen_cost,
    //  input->max_nitrogen,
    //  input->nitrogen,
    //  input->photosynthesis,
    //  input->max_photosyntetic_rate,
    //  input->lambda_k,
    //  input->optimal_stechiometric_ratio,
    //  input->phosphorus,
    //  input->min_nitrogen,
    //  input->nitrogen_uptake_sucrose_consumption,
    //  input->starch_partition_coeff,
    //  input->starch_degradation_rate);
    //
    // input->phosphorus_affinity = phosphorus_affinity(
    //     input->phosphorus_affinity,
    //     input->phosphorus_uptake,
    //     input->phosphorus_cost,
    //     input->max_phosphorus,
    //     input->phosphorus,
    //     input->photosynthesis,
    //     input->max_photosyntetic_rate,
    //     input->lambda_k,
    //     input->optimal_stechiometric_ratio,
    //     input->nitrogen,
    //     input->min_phosphorus,
    //     input->phosphorus_uptake_sucrose_consumption,
    //     input->starch_partition_coeff,
    //     input->starch_degradation_rate);

    input->transport_cost = transport_cost(
        input->sucrose_consumption_transport,
        input->respiration_frequency,
        input->sucrose,
        input->sucrose_loading_frequency,
        input->night_efficiency_starch);

    input->night_efficiency_starch = night_efficieny_starch(
        input->sucrose,
        input->max_sucrose,
        input->lambda_g,
        input->light,
        input->starch_partition_coeff);

    // 12. Rast biljke
    input->leaf_biomass += leaf_growth(
        input->lambda_sb,
        input->sucrose_root_allocation,
        input->sucrose_loading_frequency,
        input->night_efficiency_starch,
        input->leaf_biomass,
        input->leaf_deathrate,
        input->leaf_competitive_rate);

    input->root_biomass += root_growth(
        input->lambda_sb,
        input->sucrose_root_allocation,
        input->sucrose_loading_frequency,
        input->night_efficiency_starch,
        input->leaf_biomass,
        input->root_biomass,
        input->root_deathrate,
        input->root_competitive_rate);

    input->total_biomass = input->leaf_biomass + input->root_biomass;

    // 13. Dodatna kontrola saturacije

    // 14. Sigurnosna granica da ništa ne ode u NaN
    LogSimulationStep(t, input);
}

void simulate_days(int days)
{
    const int totalHours = days * 24;
    float sucroseValues[10000] = { 0 };
    float starchValues[10000] = { 0 };
    float phValues[10000] = { 0 };
    float partition[10000] = { 0 };

    Input input = {
        .light = 0.0f,
        .limitation_of_photosyntetic_rate = 0.3f,
        .max_photosyntetic_rate = 12.7f,
        .nitrogen_saturation = 12.0f,
        .phosphorus_saturation = 12.0f,
        .leaf_biomass = 2.0f, // kg/m2
        .min_leaf_biomass = 0.2f,
        .feedback_on_photosynthesis = 0.5f,
        .max_starch = 72.0f, // procenat (maks. akumulacija skroba tokom dana)
        .max_starch_degradation_rate = 6.0f, // g/m2/h
        .photoperiod = 12.0f, // sati svetla
        .optimal_stechiometric_ratio = 10.0f, // N:P
        .starch_partition_coeff = 0.5f,
        .starch_degradation_rate = 0.0f,
        .uptake_cost = 0.0f,
        .transport_cost = 0.0f,
        .respiration_frequency = 0.79f,
        .sucrose_loading_frequency = 1.98f,
        .night_efficiency_starch = 0.5f,
        .max_sucrose = 2.0f, // g/m2
        .min_sucrose = 1.3f,
        .starch_night_start = 2.0f,
        .min_starch = 0.15f,
        .lambda_sdr = 0.25f,
        .lambda_sdi = 0.10f,
        .lambda_sni = 0.08f,
        .lambda_sb = 0.00587f,
        .sucrose = 1.0f, // g/m2
        .starch = 2.0f, // g/m2
        .nitrogen = 15.3f, // %
        .phosphorus = 12.2f, // %
        .nitrogen_uptake = 0.0f,
        .phosphorus_uptake = 0.0f,
        .assimilation_cost_nitrogen = 0.08f,
        .lambda_csn = 0.0267f,
        .min_nitrogen_photosynthesis = 5.75f,
        .min_phosphorus_photosynthesis = 0.005f,
        .nutrient_conversion_parameter = 0.066f,
        .total_biomass = 3.0f, // kg/m2
        .assimilation_cost_phorphorus = 0.07f,
        .pot_nitrogen_uptake = 0.0f,
        .nitrogen_nutrient_uptake = 0.0f,
        .nitrogen_affinity = 0.5f,
        .phosphorus_affinity = 0.5f,
        .root_biomass = 1.0f, // ostatak mase
        .nitrogen_uptake_sucrose_consumption = 0.11f,
        .phosphorus_uptake_sucrose_consumption = 0.11f,
        .max_nitrogen_uptake = 6.44f,
        .max_phosphorus_uptake = 0.4f,
        .Michaelis_Menten_constant_nitrogen = 0.125f,
        .Michaelis_Menten_constant_phosphorus = 0.006736f,
        .nitrogen_soil_content = 12.0f, // mg/kg
        .lambda_k = 50.0f,
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
    float dt = 1.0f;
    int iteration = 0;
    for (int t = 0; t < totalHours; t++) {
        simulate_step(t, &input);
        sucroseValues[t] = input.sucrose;
        starchValues[t] = input.starch;
        phValues[t] = input.photosynthesis;
        partition[t] = input.total_biomass;
        printf("%d\n", t);
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
        for (int i = 1; i < totalHours; ++i) {
            DrawLine(20 + (i - 1 - scrollX) * zoom, 500 - sucroseValues[i - 1] * 20 * zoom - scrollY, 20 + (i - scrollX) * zoom, 500 - sucroseValues[i] * 20 * zoom - scrollY, RED);
            DrawLine(20 + (i - 1 - scrollX) * zoom, 500 - starchValues[i - 1] * 20 * zoom - scrollY, 20 + (i - scrollX) * zoom, 500 - starchValues[i] * 20 * zoom - scrollY, BLUE);
            DrawLine(20 + (i - 1 - scrollX) * zoom, 500 - phValues[i - 1] * 20 * zoom - scrollY, 20 + (i - scrollX) * zoom, 500 - phValues[i] * 20 * zoom - scrollY, GREEN);
            DrawLine(20 + (i - 1 - scrollX) * zoom, 500 - partition[i - 1] * 20 * zoom - scrollY, 20 + (i - scrollX) * zoom, 500 - partition[i] * 20 * zoom - scrollY, BLACK);
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

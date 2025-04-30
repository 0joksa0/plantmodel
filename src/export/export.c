#include "export/export.h"
#include <stdio.h>

void log_simulation_step(int t, Input* input)
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
                "LeafBiomass,RootBiomass,TotalBiomass,Photoperiod,"
                "LimPhotRate,MaxPhotRate,MinLeafBiomass,FeedbackOnPhot,"
                "MaxStarch,MaxStarchDegr,MinSucrose,MaxSucrose,"
                "StarchNightStart,MinStarch,LambdaSdr,LambdaSdi,LambdaSni,LambdaG,LambdaSb,"
                "LambdaCsn,AssimCostN,AssimCostP,MinNPhotosynth,MinPPhotosynth,NutrConvParam,"
                "MaxNitrogen,MinNitrogen,MaxPhosphorus,MinPhosphorus,"
                "NitSoilContent,PhosSoilContent,MaxNitUptake,MaxPhosUptake,"
                "MMConstantN,MMConstantP,RootAllocation,StochSignal,"
                "LeafDeathRate,RootDeathRate,LeafCompRate,RootCompRate\n");
        }

        fprintf(log,
            "%d,%.0f,%.4f,%.4f,"
            "%.4f,%.4f,%.4f,%.4f,%.4f,"
            "%.4f,%.4f,%.4f,%.4f,"
            "%.4f,%.4f,%.4f,"
            "%.4f,%.4f,%.4f,%.2f,"
            "%.4f,%.4f,%.4f,%.4f,"
            "%.4f,%.4f,%.4f,%.4f,"
            "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,"
            "%.4f,%.4f,%.4f,%.4f,%.4f,"
            "%.4f,%.4f,%.4f,%.4f,"
            "%.4f,%.4f,%.4f,%.4f,"
            "%.4f,%.4f,%.4f,%.4f,%.4f,"
            "%.4f,%.4f,%.4f,%.4f\n",
            t, input->light,
            input->starch_partition_coeff, input->photosynthesis,
            input->sucrose, input->starch,
            input->starch_degradation_rate,
            input->starch_partition_coeff * input->photosynthesis - input->starch_degradation_rate,
            (1.0f - input->starch_partition_coeff) * input->photosynthesis + input->starch_degradation_rate - input->uptake_cost - input->transport_cost - (input->respiration_frequency * input->sucrose) - (input->sucrose_loading_frequency * input->night_efficiency_starch),
            input->nitrogen, input->phosphorus,
            input->nitrogen_affinity, input->phosphorus_affinity,
            input->night_efficiency_starch, input->uptake_cost, input->transport_cost,
            input->leaf_biomass, input->root_biomass, input->total_biomass, input->photoperiod,
            input->limitation_of_photosyntetic_rate, input->max_photosyntetic_rate, input->min_leaf_biomass, input->feedback_on_photosynthesis,
            input->max_starch, input->max_starch_degradation_rate, input->min_sucrose, input->max_sucrose,
            input->starch_night_start, input->min_starch, input->lambda_sdr, input->lambda_sdi, input->lambda_sni, input->lambda_g, input->lambda_sb,
            input->lambda_csn, input->assimilation_cost_nitrogen, input->assimilation_cost_phorphorus, input->min_nitrogen_photosynthesis, input->min_phosphorus_photosynthesis, input->nutrient_conversion_parameter,
            input->max_nitrogen, input->min_nitrogen, input->max_phosphorus, input->min_phosphorus,
            input->nitrogen_soil_content, input->phosphorus_soil_content, input->max_nitrogen_uptake, input->max_phosphorus_uptake,
            input->Michaelis_Menten_constant_nitrogen, input->Michaelis_Menten_constant_phosphorus, input->sucrose_root_allocation, input->stochiometric_signal,
            input->leaf_deathrate, input->root_deathrate, input->leaf_competitive_rate, input->root_competitive_rate);

        fclose(log);
    }
}


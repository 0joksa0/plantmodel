#include "export/export.h"
#include "model/model.h"
#include <stdio.h>

void log_simulation_step(int days,real_t t, Input* input)
{
    /* TODO: Keep file handle open for the full simulation; reopening on every step is expensive and can lose data on interruption. */

    char filename[64];
    snprintf(filename, sizeof(filename), "sim_log_%.0fh_%dd.csv", (double)input->photo.photoperiod, days);
    FILE* log = fopen(filename, t == 0 ? "w" : "a");
    if (log) {
        if (t == REAL(0.0)) {
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
            "%.4f,%.0f,%.4f,%.4f,"
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
            t, input->photo.light,
            input->carbohydrates.starch_partition_coeff, input->photo.photosynthesis,
            input->carbohydrates.sucrose, input->carbohydrates.starch,
            input->carbohydrates.starch_degradation_rate,
            starch_production(input->photo.photosynthesis, input->carbohydrates.starch_partition_coeff, input->carbohydrates.starch_degradation_rate),
            sucrose_production(input->carbohydrates.starch_partition_coeff, input->photo.photosynthesis, input->carbohydrates.starch_degradation_rate, input->carbohydrates.uptake_cost, input->carbohydrates.transport_cost, input->carbohydrates.respiration_frequency, input->carbohydrates.sucrose, input->carbohydrates.sucrose_loading_frequency, input->carbohydrates.night_efficiency_starch),
            input->nutrients.nitrogen, input->nutrients.phosphorus,
            input->nutrients.nitrogen_affinity, input->nutrients.phosphorus_affinity,
            input->carbohydrates.night_efficiency_starch, input->carbohydrates.uptake_cost, input->carbohydrates.transport_cost,
            input->growth.leaf_biomass, input->growth.root_biomass, input->growth.total_biomass, input->photo.photoperiod,
            input->photo.limitation_of_photosynthetic_rate, input->photo.max_photosynthetic_rate, input->photo.min_leaf_biomass, input->photo.feedback_on_photosynthesis,
            input->carbohydrates.max_starch, input->carbohydrates.max_starch_degradation_rate, input->carbohydrates.min_sucrose, input->carbohydrates.max_sucrose,
            input->carbohydrates.starch_night_start, input->carbohydrates.min_starch, input->carbohydrates.lambda_sdr, input->carbohydrates.lambda_sdi, input->carbohydrates.lambda_sni, input->carbohydrates.lambda_g, input->growth.lambda_sb,
            input->nutrients.lambda_csn, input->nutrients.assimilation_cost_nitrogen, input->nutrients.assimilation_cost_phosphorus, input->nutrients.min_nitrogen_photosynthesis, input->nutrients.min_phosphorus_photosynthesis, input->nutrients.nutrient_conversion_parameter,
            input->nutrients.max_nitrogen, input->nutrients.min_nitrogen, input->nutrients.max_phosphorus, input->nutrients.min_phosphorus,
            input->nutrients.nitrogen_soil_content, input->nutrients.phosphorus_soil_content, input->nutrients.max_nitrogen_uptake, input->nutrients.max_phosphorus_uptake,
            input->nutrients.Michaelis_Menten_constant_nitrogen, input->nutrients.Michaelis_Menten_constant_phosphorus, input->growth.sucrose_root_allocation, input->nutrients.stoichiometric_signal,
            input->growth.leaf_deathrate, input->growth.root_deathrate, input->growth.leaf_competitive_rate, input->growth.root_competitive_rate);

        fclose(log);
    }
}

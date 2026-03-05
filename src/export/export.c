#include "export/export.h"
#include "model/model.h"
#include <stdio.h>

void log_simulation_step(int days,real_t t, Input* input)
{
    /* TODO: Keep file handle open for the full simulation; reopening on every step is expensive and can lose data on interruption. */

    char filename[64];
    snprintf(filename, sizeof(filename), "sim_log_%.0fh_%dd.csv", (double)input->core.photoperiod, days);
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
            t, input->core.light,
            input->core.starch_partition_coeff, input->core.photosynthesis,
            input->core.sucrose, input->core.starch,
            input->core.starch_degradation_rate,
            starch_production(input->core.photosynthesis, input->core.starch_partition_coeff, input->core.starch_degradation_rate),
            sucrose_production(input->core.starch_partition_coeff, input->core.photosynthesis, input->core.starch_degradation_rate, input->core.uptake_cost, input->core.transport_cost, input->core.respiration_frequency, input->core.sucrose, input->core.sucrose_loading_frequency, input->core.night_efficiency_starch),
            input->core.nitrogen, input->core.phosphorus,
            input->core.nitrogen_affinity, input->core.phosphorus_affinity,
            input->core.night_efficiency_starch, input->core.uptake_cost, input->core.transport_cost,
            input->core.leaf_biomass, input->core.root_biomass, input->core.total_biomass, input->core.photoperiod,
            input->core.limitation_of_photosynthetic_rate, input->core.max_photosynthetic_rate, input->core.min_leaf_biomass, input->core.feedback_on_photosynthesis,
            input->core.max_starch, input->core.max_starch_degradation_rate, input->core.min_sucrose, input->core.max_sucrose,
            input->core.starch_night_start, input->core.min_starch, input->core.lambda_sdr, input->core.lambda_sdi, input->core.lambda_sni, input->core.lambda_g, input->core.lambda_sb,
            input->core.lambda_csn, input->core.assimilation_cost_nitrogen, input->core.assimilation_cost_phosphorus, input->core.min_nitrogen_photosynthesis, input->core.min_phosphorus_photosynthesis, input->core.nutrient_conversion_parameter,
            input->core.max_nitrogen, input->core.min_nitrogen, input->core.max_phosphorus, input->core.min_phosphorus,
            input->core.nitrogen_soil_content, input->core.phosphorus_soil_content, input->core.max_nitrogen_uptake, input->core.max_phosphorus_uptake,
            input->core.Michaelis_Menten_constant_nitrogen, input->core.Michaelis_Menten_constant_phosphorus, input->core.sucrose_root_allocation, input->core.stoichiometric_signal,
            input->core.leaf_deathrate, input->core.root_deathrate, input->core.leaf_competitive_rate, input->core.root_competitive_rate);

        fclose(log);
    }
}

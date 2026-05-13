#include "model/model.h"

real_t leaf_growth(
    real_t lambda_sb,
    real_t sucrose_root_allocation,
    real_t sucrose_loading_frequency,
    real_t night_efficiency_starch,
    real_t leaf_biomass,
    real_t leaf_death_rate,
    real_t leaf_competitive_rate)
{
    return (lambda_sb * (REAL(1.0) - sucrose_root_allocation) * sucrose_loading_frequency * night_efficiency_starch * leaf_biomass) - (leaf_death_rate * leaf_biomass) - (leaf_competitive_rate * leaf_biomass * leaf_biomass);
}

real_t leaf_growth_f(
    real_t t,
    real_t leaf_biomass,
    void* params)
{
    Input* input = (Input*)params;
    (void)t;

    return leaf_growth(
        input->growth.lambda_sb,
        input->growth.sucrose_root_allocation,
        input->carbohydrates.sucrose_loading_frequency,
        input->carbohydrates.night_efficiency_starch,
        leaf_biomass,
        input->growth.leaf_deathrate,
        input->growth.leaf_competitive_rate);
}

real_t root_growth(
    real_t lambda_sb,
    real_t sucrose_root_allocation,
    real_t sucrose_loading_frequency,
    real_t night_efficiency_starch,
    real_t leaf_biomass,
    real_t root_biomass,
    real_t root_death_rate,
    real_t root_competitive_rate)
{
    return (lambda_sb * sucrose_root_allocation * sucrose_loading_frequency * night_efficiency_starch * leaf_biomass) - (root_death_rate * root_biomass) - (root_competitive_rate * root_biomass * root_biomass);
}

real_t root_growth_f(
    real_t t,
    real_t root_biomass,
    void* params)
{

    Input* input = (Input*)params;
    (void)t;
    return root_growth(
        input->growth.lambda_sb,
        input->growth.sucrose_root_allocation,
        input->carbohydrates.sucrose_loading_frequency,
        input->carbohydrates.night_efficiency_starch,
        input->growth.leaf_biomass,
        root_biomass,
        input->growth.root_deathrate,
        input->growth.root_competitive_rate);
}

real_t sucrose_root_allocation(
    real_t sucrose_root_allocation,
    real_t nitrogen_affinity,
    real_t phosphorus_affinity,
    real_t stoichiometric_signal,
    real_t sucrose,
    real_t min_sucrose,
    real_t nitrogen,
    real_t min_nitrogen,
    real_t phosphorus,
    real_t min_phosphorus)
{
    real_t first = (REAL(1.0) - sucrose_root_allocation) * ((nitrogen_affinity * stoichiometric_signal) + ((REAL(1.0) - stoichiometric_signal) * phosphorus_affinity));
    real_t second = sucrose_root_allocation * (((nitrogen * stoichiometric_signal) / (nitrogen + min_nitrogen)) + ((phosphorus * (REAL(1.0) - stoichiometric_signal)) / (phosphorus + min_phosphorus)) + (min_sucrose / (min_sucrose + sucrose)));
    return first - second;
}

real_t sucrose_root_allocation_f(
    real_t t,
    real_t sucrose_root_allocation_p,
    void* params)
{
    Input* input = (Input*)params;
    (void)t;
    return sucrose_root_allocation(
        sucrose_root_allocation_p,
        input->nutrients.nitrogen_affinity,
        input->nutrients.phosphorus_affinity,
        input->nutrients.stoichiometric_signal,
        input->carbohydrates.sucrose,
        input->carbohydrates.min_sucrose,
        input->nutrients.nitrogen,
        input->nutrients.min_nitrogen,
        input->nutrients.phosphorus,
        input->nutrients.min_phosphorus);
}

real_t lambda_sb_f(real_t lambda_sb, real_t nitrogen_soil, real_t phosphorus_soil)
{
    real_t n = nitrogen_soil / REAL(12.0);
    real_t p = phosphorus_soil / REAL(0.15);
    real_t g = -(REAL(0.5) * (RPOW(n, 2) + RPOW(p, 2))) + (n + p);
    real_t R_minus = RMAX(g, 0);
    if (RMIN(n, p) >= 1) {
        return lambda_sb;
    } else {
        return lambda_sb * R_minus;
    }
}

real_t R_plus(real_t nitrogen_soil, real_t phosphorus_soil)
{
    real_t n = nitrogen_soil / REAL(12.0);
    real_t p = phosphorus_soil / REAL(0.15);

    real_t R_plus = (REAL(69.0) * (RPOW(n, 2) + RPOW(p, 2))) - REAL(138.0) * (n + p) + REAL(139);
    if (RMIN(n, p) >= 1) {
        return 1;
    }
    return R_plus;
}

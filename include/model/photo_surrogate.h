#ifndef PHOTO_SURROGATE_H
#define PHOTO_SURROGATE_H

double photo_surrogate_should_use(double light_PAR);

double photo_surrogate_predict_an(
    double light_PAR,
    double leaf_biomass,
    double stomatal_conductance);

double photo_surrogate_predict_if_applicable(
    double light_PAR,
    double leaf_biomass,
    double stomatal_conductance,
    int *used_surrogate);

#endif

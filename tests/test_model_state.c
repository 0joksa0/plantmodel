#include <criterion/criterion.h>
#include <solver.h>

#include "model/algebraic.h"
#include "model/input.h"
#include "model/environment.h"
#include "model/outputs.h"
#include "model/parameters.h"
#include "model/state.h"
#include "simulation/dynamics.h"

Test(model_state, pack_unpack_round_trip_preserves_core_state)
{
    Input input = generate_input();
    real_t state[X_DIM] = {0};
    Input round_trip = input;

    simulation_pack_state(state, &input);
    round_trip.carbohydrates.starch = REAL(-1.0);
    round_trip.carbohydrates.sucrose = REAL(-1.0);
    round_trip.growth.total_biomass = REAL(-1.0);
    simulation_unpack_state(&round_trip, state);

    cr_assert(round_trip.carbohydrates.starch == input.carbohydrates.starch);
    cr_assert(round_trip.carbohydrates.sucrose == input.carbohydrates.sucrose);
    cr_assert(round_trip.growth.total_biomass == input.growth.total_biomass);
}

Test(model_state, algebraic_compute_matches_reference_fields)
{
    Input input = generate_input();

    simulation_compute_algebraic(&input, REAL(6.0), input.carbohydrates.starch);

    cr_assert(input.photo.photosynthesis > REAL(0.0));
    cr_assert(input.growth.total_biomass > REAL(0.0));
    cr_assert(input.photo.nitrogen_saturation >= REAL(0.0));
    cr_assert(input.photo.phosphorus_saturation >= REAL(0.0));
}

Test(model_state, typed_state_round_trip_matches_solver_vector_layout)
{
    Input input = generate_input();
    PlantState typed = plant_state_from_input(&input);
    real_t state[X_DIM] = {0};
    PlantState unpacked = {0};

    plant_state_pack(state, &typed);
    plant_state_unpack(&unpacked, state);

    cr_assert(unpacked.starch == typed.starch);
    cr_assert(unpacked.sucrose == typed.sucrose);
    cr_assert(unpacked.root_biomass == typed.root_biomass);
}

Test(model_state, typed_outputs_match_legacy_algebraic_path)
{
    Input input = generate_input();
    PlantState state = plant_state_from_input(&input);
    PlantParameters params = plant_parameters_from_input(&input);
    PlantEnvironment env = plant_environment_from_input(&input, REAL(6.0));
    PlantOutputs outputs = {0};

    plant_compute_outputs(&state, &params, &env, input.carbohydrates.starch, &outputs);
    simulation_compute_algebraic(&input, REAL(6.0), input.carbohydrates.starch);

    cr_assert(outputs.photosynthesis == input.photo.photosynthesis);
    cr_assert(outputs.total_biomass == input.growth.total_biomass);
    cr_assert(outputs.nitrogen_saturation == input.photo.nitrogen_saturation);
    cr_assert(outputs.phosphorus_saturation == input.photo.phosphorus_saturation);
}

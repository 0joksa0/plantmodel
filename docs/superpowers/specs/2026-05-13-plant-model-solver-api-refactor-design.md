# Plant Model Solver API Refactor Design

**Date:** 2026-05-13

## Goal

Refactor the plant model so it becomes a first-class client of the new solver/framework API, with explicit separation between state, parameters, environment inputs, outputs, runtime context, and simulation orchestration, while preserving the existing mathematics and numerical behavior.

## Constraints

- Do not change plant mathematics or intended biological behavior.
- Preserve current solver-backed simulation capabilities: headless runs, GUI, CSV export, observers, runtime stats.
- Keep `Input` temporarily as a legacy configuration and initialization boundary.
- Prioritize work aligned with `#shared-v02` sprint tasks:
  - state vector extraction
  - generic model interface
  - verifiable solver behavior

## Current Problems

- The code formally uses `ModelInterface` and `model_run`, but the internal execution path is still centered on mutable `Input`.
- `ModelDescriptor` exposes only state fields; parameters, inputs, and outputs are not modeled.
- RHS and algebraic calculations reconstruct and mutate `Input` snapshots instead of operating on explicit model buffers.
- Model-specific runtime state (`starch_night_start`, `was_light`) is hidden inside simulation/observer glue rather than owned by the model runtime context.
- `simulation` still knows too much about plant state layout and model internals.

## Target Architecture

The refactor will follow a strangler approach:

1. Introduce explicit model data types without changing behavior.
2. Move state layout and mathematical domains into focused model modules.
3. Promote the plant model to a complete solver-framework model with state, parameter, input, and output descriptors.
4. Reduce `simulation` to orchestration, sampling, GUI, CSV, and lifecycle concerns.

Execution path:

`config/init -> state/parameter/environment setup -> model outputs/algebraic compute -> RHS compute -> solver runtime -> snapshot-based observers/export`

## Data Model

### `PlantState`

Owns only dynamic solver state. This becomes the authoritative mapping for the solver state vector.

Expected fields:

- `starch_partition_coeff`
- `starch`
- `sucrose`
- `nitrogen_affinity`
- `phosphorus_affinity`
- `nitrogen`
- `phosphorus`
- `sucrose_root_allocation`
- `leaf_biomass`
- `root_biomass`

### `PlantParameters`

Owns model constants and slowly changing coefficients currently embedded in `Input.core`. This includes growth, uptake, stoichiometry, respiration, and gas-exchange configuration values that are not dynamic solver state.

### `PlantEnvironment`

Owns exogenous or time-dependent model inputs.

Expected categories:

- time-of-day derived light state
- PAR / light-related environment values
- gas-exchange environment values
- soil nutrient environment values
- photoperiod and schedule-derived environment data

The first iteration may still populate this through a legacy adapter from `Input`.

### `PlantOutputs`

Owns algebraic and diagnostic values derived from state + parameters + environment.

Expected categories:

- photosynthesis
- nutrient saturations
- starch degradation rate
- uptake and transport costs
- nutrient uptake intermediates
- total biomass
- any other non-state values currently computed inside `simulation_compute_algebraic`

### `PlantSnapshot`

Read-only composition of:

- `PlantState`
- `PlantParameters`
- `PlantEnvironment`
- `PlantOutputs`

This becomes the common view used by observers, CSV export, GUI sampling, and tests.

### `PlantModelCtx`

Minimal runtime context owned by the model execution path.

Expected contents:

- reference to parameters
- reference or callback for environment update
- runtime phase-transition cache:
  - `starch_night_start`
  - `was_light`

This is model-specific runtime state, not a solver concern.

## Module Layout

### New or repurposed core model files

- `include/model/state.h`, `src/model/state.c`
  - state indices
  - state vector pack/unpack
  - conversions between typed state and solver buffers

- `include/model/parameters.h`, `src/model/parameters.c`
  - parameter type
  - legacy adapter from `Input`

- `include/model/environment.h`, `src/model/environment.c`
  - environment type
  - light/day schedule update
  - legacy adapter hooks

- `include/model/outputs.h`, `src/model/outputs.c`
  - outputs type
  - typed helpers for populating outputs from current model inputs

- `include/model/context.h`
  - `PlantModelCtx`

- `include/model/snapshot.h`
  - `PlantSnapshot`

### Domain math modules

Mathematics is preserved; functions are reorganized by domain only.

- `include/model/photosynthesis.h`, `src/model/photosynthesis.c`
- `include/model/carbon.h`, `src/model/carbon.c`
- `include/model/nutrients.h`, `src/model/nutrients.c`
- `include/model/growth.h`, `src/model/growth.c`
- `include/model/algebraic.h`, `src/model/algebraic.c`
- `include/model/environment.h`, `src/model/environment.c`

### Public model API

- `include/model/model.h`, `src/model/model.c`
  - expose `ModelDescriptor`
  - expose `ModelInterface`
  - implement output computation entrypoint
  - implement RHS entrypoint

### Simulation layer

- `src/simulation/`
  - owns run configuration and orchestration only
  - owns observer chain, GUI thread, CSV lifecycle, logging, sampling
  - no direct ownership of plant mathematics or state indexing

## Runtime Flow

1. `config_load_file(...)` fills legacy `Input` and `SimulationConfig`.
2. Adapter layer builds:
   - initial `PlantState`
   - `PlantParameters`
   - initial `PlantEnvironment`
3. `simulation` initializes `ModelBuffers` with state, parameters, inputs, and outputs buffers.
4. Before each solver evaluation:
   - environment is refreshed for time `t`
   - model runtime phase-transition helper updates transition state if needed
   - outputs are computed
   - RHS is computed from state + parameters + environment + outputs
5. Observers and exporters consume a `PlantSnapshot`, not reconstructed `Input`.
6. At run end, the final solver state is unpacked back into legacy-facing structures as needed.

## Night Transition Handling

The current `night_start` logic must become explicit and centralized.

Introduce a model-owned runtime helper, for example:

- `plant_runtime_update_phase_transition(...)`

Responsibilities:

- detect day-to-night transition from current environment/light state
- capture `starch_night_start` from the current solver state
- update `was_light`

This helper belongs to the model runtime context and is invoked from one controlled execution path. It must not remain split across simulation, RHS setup, and observers.

## Descriptor Requirements

The final `ModelDescriptor` must expose:

- state fields
- parameter fields
- input fields (`PlantEnvironment`)
- output fields (`PlantOutputs`)

The current state-only descriptor is an intermediate state and must be replaced.

## Migration Strategy

### Phase 1: Introduce explicit model types

Add `PlantState`, `PlantParameters`, `PlantEnvironment`, `PlantOutputs`, `PlantSnapshot`, and `PlantModelCtx` without changing behavior.

### Phase 2: Extract state vector ownership

Move state indexing and pack/unpack logic out of `simulation` into the model state layer.

### Phase 3: Split mathematical domains

Reorganize existing mathematical functions into domain modules without changing formulas.

### Phase 4: Extract outputs/algebraic layer

Replace `simulation_compute_algebraic(...)` with a model-owned outputs computation layer.

### Phase 5: Introduce environment provider path

Route light and other exogenous inputs through `PlantEnvironment`, using legacy adapters where needed.

### Phase 6: Promote full model descriptor

Expose states, parameters, inputs, and outputs to the solver framework.

### Phase 7: Convert simulation consumers to snapshots

Switch observers, CSV export, GUI sampling, and similar code to `PlantSnapshot`.

### Phase 8: Shrink legacy `Input`

Keep `Input` only as config/init boundary, then progressively eliminate deeper dependencies once the new path is stable.

## Verification Strategy

Refactor verification is regression-oriented.

### Existing gates to preserve

- current photoperiod/RGR tests remain active
- current end-to-end solver-backed simulation path remains active

### New regression gates

- state pack/unpack equivalence tests
- algebraic/output snapshot equivalence tests
- RHS equivalence tests for fixed state/time inputs
- observer/sampling snapshot tests for key exported values

The primary evidence goal is architectural equivalence, not new biological behavior.

## Non-Goals

- no mathematical retuning
- no deliberate formula changes
- no new biological features
- no full removal of legacy `Input` in the first implementation pass

## Risks

- accidental math drift during function relocation
- state mapping bugs during typed-state extraction
- hidden dependencies on mutable `Input` fields in observer/export code
- incorrect handling of day/night transition cache when moving logic into model context

## Success Criteria

- plant model executes through a complete solver-framework model interface
- `simulation` no longer owns plant state indexing or algebraic model logic
- observers/export/GUI read from snapshots rather than reconstructed `Input`
- regression tests show no material change in mathematical behavior

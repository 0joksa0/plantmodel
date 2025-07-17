
# PlantModel – Digital Twin Plant Growth Simulator in C

## Overview

**PlantModel** is a modular simulation engine written in C designed to model the physiological growth of a plant using a digital twin paradigm. The simulator integrates environmental and internal biological parameters to simulate plant development over time (e.g., leaf/stem/root growth, nitrogen and phosphorus use).

This engine relies on [solver](https://github.com/0joksa0/solver), a custom high-precision numerical ODE solver also written in C. Combined with real-time Raylib rendering and structured data export, PlantModel is suitable for scientific simulations, digital agriculture, and academic research in crop modeling.

It was developed as part of a master's thesis focused on plant digital twins and simulation-based decision support.

## Features

- Numerical simulation of plant physiology over time
- Raylib-powered GUI with real-time biomass visualization
- Modular structure: input model, solver integration, export
- Compatible with scientific data formats (CSV)
- Designed for easy parameter modification and future model expansion

## Project Structure

```
plantModel/
├── src/
│   └── main.c
├── gui/
│   └── plot.c
├── export/
│   └── export.c
├── simulation/
│   └── simulation.c
├── model/
│   └── input.c
│   └── model.c
├── include/
├── CMakeLists.txt
```
```

## Build Instructions

### Requirements

- [`solver`](https://github.com/0joksa0/solver) static library (`solver_double`)
- `raylib` (with `pkg-config`)
- System libs: `libm`, `pthread`, `libdl`, `GL`

### Build

```bash
cmake -B build 
cmake --build build
```

Binary will be located at: `build/simulator`

## Usage

```bash
./simulator
```

The program:
- Initializes plant inputs (light, NDVI, nutrients)
- Runs a 29-day simulation (`simulate_days`)
- Displays a biomass-over-time graph
- Optionally saves results to CSV

## Future Work

- Integrate environmental inputs (e.g. NDVI from satellite)
- Add root modeling and water stress response
- Support for plant-type parameterization
- Real-time control/feedback via sensors
- Calibration from experimental data



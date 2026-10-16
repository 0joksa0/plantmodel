cmake_minimum_required(VERSION 3.20)
project(digital_twin_c C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ---------------- kompilacija ----------------
# TODO: Split compile flags by build type (`Debug` vs `Release`) and avoid forcing `-ffast-math` globally.
 add_compile_options(-Wall -Wextra -O3 -ffast-math)

find_package(OpenMP)

# ---------------- fajlovi --------------------
# Feature definitions for the revised model.  None of these are optional -- each
# one was adopted on a measurement, and dropping any of them changes the reported
# results:
#   SINK_ACTIVATION            sugar-induced sink term; without it RGR MAPE goes
#                              from 8.6% to 36.7%
#   HAS_LIGHT_PROFILE          square-wave chamber light; without it --chamber
#   HAS_DAYNIGHT_TEMP          silently runs a half-sine at a single temperature
#   STARCH_FEEDBACK_FARQUHAR   routes the starch limitation into the Farquhar
#                              path (it is inert there at lambda_c = 1, which is
#                              why lambda_c is not a fitted parameter)
#   RESPIRATION_ALL_HOURS      maintenance respiration runs day and night, per
#                              Tedone Eq. 2.2.30, rather than in darkness only
# They lived in CMAKE_C_FLAGS in the development sandbox, so a plain configure of
# this tree reproduced 9.2% instead of 8.6% until they were made explicit here.
add_compile_definitions(
    PLANTMODEL_HAS_LIGHT_PROFILE
    PLANTMODEL_HAS_DAYNIGHT_TEMP
    PLANTMODEL_SINK_ACTIVATION
    PLANTMODEL_STARCH_FEEDBACK_FARQUHAR
    PLANTMODEL_RESPIRATION_ALL_HOURS
)
add_compile_definitions(SOLVER_REAL_AS=double )
add_compile_definitions(SOLVER_REAL_IS_DOUBLE)
file(GLOB_RECURSE SRC CONFIGURE_DEPENDS src/*.c)
add_executable(simulator ${SRC})
target_include_directories(simulator PRIVATE ${PROJECT_SOURCE_DIR}/include)

# ---------------- solver_double --------------
# 1) nađi statičku ili deljenu biblioteku
set(SOLVER_PREFIX "/usr/local" CACHE PATH "Prefix where solver headers and libraries are installed")
find_library(SOLVER_DOUBLE_PATH
             NAMES solver_double
             PATHS "${SOLVER_PREFIX}/lib"
             NO_DEFAULT_PATH)

if(NOT SOLVER_DOUBLE_PATH)
    find_library(SOLVER_DOUBLE_PATH
                 NAMES solver_double
                 PATHS /usr/local/lib /usr/lib
                 REQUIRED)
endif()

find_path(SOLVER_INCLUDE_PATH
          NAMES solver.h
          PATHS "${SOLVER_PREFIX}/include"
          NO_DEFAULT_PATH)

if(NOT SOLVER_INCLUDE_PATH)
    find_path(SOLVER_INCLUDE_PATH
              NAMES solver.h
              PATHS /usr/local/include /usr/include
              REQUIRED)
endif()

# 2) napravi IMPORTED target
add_library(solver_double STATIC IMPORTED)
set_target_properties(solver_double PROPERTIES
    IMPORTED_LOCATION             ${SOLVER_DOUBLE_PATH}
    INTERFACE_INCLUDE_DIRECTORIES "${SOLVER_INCLUDE_PATH}"
)

# ---------------- raylib ----------------------
find_package(PkgConfig REQUIRED)
pkg_check_modules(RAYLIB REQUIRED IMPORTED_TARGET raylib)

# ---------------- sistemske libs --------------
find_library(LIBM   m       REQUIRED)
find_library(LIBDL  dl      REQUIRED)
find_library(LIBPTH pthread REQUIRED)
find_library(LIBGL  GL      REQUIRED)

# ---------------- link ------------------------

target_link_libraries(simulator
    PRIVATE
        solver_double 
        PkgConfig::RAYLIB
        ${LIBM}
        ${LIBPTH}
        ${LIBDL}
        ${LIBGL}
)

list(FILTER SRC EXCLUDE REGEX ".*/main\\.c$")
add_executable(sweep_runner ${SRC} tools/sweep_runner.c)
target_include_directories(sweep_runner PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(sweep_runner
    PRIVATE
        solver_double
        PkgConfig::RAYLIB
        ${LIBM}
        ${LIBPTH}
        ${LIBDL}
        ${LIBGL}
)

add_executable(validation_runner ${SRC} tools/validation_runner.c)

add_executable(model_report ${SRC} tools/model_report.c)
target_include_directories(model_report PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(model_report PRIVATE solver_double PkgConfig::RAYLIB ${LIBM} ${LIBPTH} ${LIBDL} ${LIBGL})

# Diel trajectory export for the carbohydrate panels.  Takes --shared-lc and
# --shared-lsni so it reports the same parameterisation as model_report; without
# them it falls back to the per-photoperiod literature values, which is how the
# carbohydrate figures once described a different model from the RGR table.
add_executable(diel_trace ${SRC} tools/diel_trace.c)
target_include_directories(diel_trace PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(diel_trace PRIVATE solver_double PkgConfig::RAYLIB ${LIBM} ${LIBPTH} ${LIBDL} ${LIBGL})

add_executable(nsweep ${SRC} tools/nsweep.c)
target_include_directories(nsweep PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(nsweep PRIVATE solver_double PkgConfig::RAYLIB ${LIBM} ${LIBPTH} ${LIBDL} ${LIBGL})
target_include_directories(validation_runner PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(validation_runner
    PRIVATE
        solver_double
        PkgConfig::RAYLIB
        ${LIBM}
        ${LIBPTH}
        ${LIBDL}
        ${LIBGL}
)

add_executable(recalibration ${SRC} tools/recalibration.c)
target_include_directories(recalibration PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_compile_options(recalibration PRIVATE -fopenmp)
target_link_libraries(recalibration
    PRIVATE
        solver_double
        PkgConfig::RAYLIB
        ${LIBM}
        ${LIBPTH}
        ${LIBDL}
        ${LIBGL}
        -fopenmp
)

add_executable(sensitivity_runner ${SRC} tools/sensitivity_runner.c)
target_include_directories(sensitivity_runner PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(sensitivity_runner
    PRIVATE
        solver_double
        PkgConfig::RAYLIB
        ${LIBM}
        ${LIBPTH}
        ${LIBDL}
        ${LIBGL}
)

add_executable(ablation_runner ${SRC} tools/ablation_runner.c)

add_executable(lg_sweep ${SRC} tools/lg_sweep.c)
target_include_directories(lg_sweep PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(lg_sweep PRIVATE solver_double PkgConfig::RAYLIB ${LIBM} ${LIBPTH} ${LIBDL} ${LIBGL})

add_executable(fj_diel_runner ${SRC} tools/fj_diel_runner.c)
target_include_directories(fj_diel_runner PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(fj_diel_runner PRIVATE solver_double PkgConfig::RAYLIB ${LIBM} ${LIBPTH} ${LIBDL} ${LIBGL})

add_executable(vacuole_check ${SRC} tools/vacuole_check.c)
target_include_directories(vacuole_check PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(vacuole_check PRIVATE solver_double PkgConfig::RAYLIB ${LIBM} ${LIBPTH} ${LIBDL} ${LIBGL})

add_executable(fj_rgr ${SRC} tools/fj_rgr.c)
target_include_directories(fj_rgr PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(fj_rgr PRIVATE solver_double PkgConfig::RAYLIB ${LIBM} ${LIBPTH} ${LIBDL} ${LIBGL})

add_executable(solver_verification ${SRC} tools/solver_verification.c)
target_include_directories(solver_verification PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(solver_verification PRIVATE solver_double PkgConfig::RAYLIB ${LIBM} ${LIBPTH} ${LIBDL} ${LIBGL})

add_executable(crossval_runner ${SRC} tools/crossval_runner.c)
target_include_directories(crossval_runner PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(crossval_runner PRIVATE solver_double PkgConfig::RAYLIB ${LIBM} ${LIBPTH} ${LIBDL} ${LIBGL})

add_executable(horizon_study ${SRC} tools/horizon_study.c)
target_include_directories(horizon_study PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(horizon_study PRIVATE solver_double PkgConfig::RAYLIB ${LIBM} ${LIBPTH} ${LIBDL} ${LIBGL})

add_executable(carbon_balance ${SRC} tools/carbon_balance.c)
target_include_directories(carbon_balance PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(carbon_balance PRIVATE solver_double PkgConfig::RAYLIB ${LIBM} ${LIBPTH} ${LIBDL} ${LIBGL})

add_executable(carbon_audit ${SRC} tools/carbon_audit.c)
target_include_directories(carbon_audit PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(carbon_audit PRIVATE solver_double PkgConfig::RAYLIB ${LIBM} ${LIBPTH} ${LIBDL} ${LIBGL})

add_executable(convergence_study ${SRC} tools/convergence_study.c)
target_include_directories(convergence_study PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(convergence_study PRIVATE solver_double PkgConfig::RAYLIB ${LIBM} ${LIBPTH} ${LIBDL} ${LIBGL})

add_executable(fj_sucrose_calibrate ${SRC} tools/fj_sucrose_calibrate.c)
target_include_directories(fj_sucrose_calibrate PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(fj_sucrose_calibrate PRIVATE solver_double PkgConfig::RAYLIB ${LIBM} ${LIBPTH} ${LIBDL} ${LIBGL})

add_executable(env_trajectory ${SRC} tools/env_trajectory.c)
target_include_directories(env_trajectory PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(env_trajectory PRIVATE solver_double PkgConfig::RAYLIB ${LIBM} ${LIBPTH} ${LIBDL} ${LIBGL})

add_executable(env_sweep ${SRC} tools/env_sweep.c)
target_include_directories(env_sweep PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(env_sweep PRIVATE solver_double PkgConfig::RAYLIB ${LIBM} ${LIBPTH} ${LIBDL} ${LIBGL})

add_executable(oat_sweep ${SRC} tools/oat_sweep.c)
target_include_directories(oat_sweep PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(oat_sweep PRIVATE solver_double PkgConfig::RAYLIB ${LIBM} ${LIBPTH} ${LIBDL} ${LIBGL})

add_executable(multi_2d_sweep ${SRC} tools/multi_2d_sweep.c)
target_include_directories(multi_2d_sweep PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(multi_2d_sweep PRIVATE solver_double PkgConfig::RAYLIB ${LIBM} ${LIBPTH} ${LIBDL} ${LIBGL})
target_include_directories(ablation_runner PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(ablation_runner
    PRIVATE
        solver_double
        PkgConfig::RAYLIB
        ${LIBM}
        ${LIBPTH}
        ${LIBDL}
        ${LIBGL}
)
if(OpenMP_C_FOUND)
    target_link_libraries(ablation_runner PRIVATE OpenMP::OpenMP_C)
endif()

message(STATUS "Solver at : ${SOLVER_DOUBLE_PATH}")
message(STATUS "Sources   : ${SRC}")


# ---------------- Criterion test framework ----------------
pkg_check_modules(CRITERION REQUIRED IMPORTED_TARGET criterion)

# Only the top-level Criterion suites; the standalone probes in
# tests/numerical/ each carry their own main() and are registered separately.
file(GLOB TEST_SRC CONFIGURE_DEPENDS tests/*.c)

if(TEST_SRC)
    add_executable(unit_tests ${SRC} ${TEST_SRC})
    target_include_directories(unit_tests PRIVATE ${PROJECT_SOURCE_DIR}/include)
    target_link_libraries(unit_tests
        PRIVATE
            PkgConfig::CRITERION
            solver_double
            PkgConfig::RAYLIB
            ${LIBM}
            ${LIBPTH}
            ${LIBDL}
            ${LIBGL}
    )
    enable_testing()
    add_test(NAME all_tests COMMAND unit_tests --jobs 1 --verbose=1 --timeout 300)
    set_tests_properties(all_tests PROPERTIES TIMEOUT 360)
endif()


# ---------------- standalone numerical probes ----------------
# Each file in tests/numerical/ is a self-contained program with its own main()
# and a gui_main_thread() stub, so it links the model sources WITHOUT
# src/gui/plot.c (raylib GUI). They resolve config/ relative to the repo root.
enable_testing()

# model/simulation sources minus the GUI (plot.c defines the real gui_main_thread)
set(SRC_NOGUI ${SRC})
list(FILTER SRC_NOGUI EXCLUDE REGEX ".*/gui/.*")

# solver-only probe: no model sources needed
add_executable(numprobe_integrator_consistency tests/numerical/integrator_consistency.c)
target_include_directories(numprobe_integrator_consistency PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(numprobe_integrator_consistency PRIVATE solver_double ${LIBM})
add_test(NAME integrator_consistency COMMAND numprobe_integrator_consistency)

# freestanding probe: no model sources, no solver
add_executable(numprobe_nan_guard tests/numerical/nan_guard_under_fastmath.c)
target_link_libraries(numprobe_nan_guard PRIVATE ${LIBM})
add_test(NAME nan_guard_under_fastmath COMMAND numprobe_nan_guard)

# full-model probes: link SRC_NOGUI, run from repo root for config/ access
add_executable(numprobe_integrator_stall ${SRC_NOGUI} tests/numerical/integrator_stall_threshold.c)
target_include_directories(numprobe_integrator_stall PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(numprobe_integrator_stall PRIVATE solver_double ${LIBM} ${LIBPTH} ${LIBDL})
add_test(NAME integrator_stall_threshold COMMAND numprobe_integrator_stall)
set_tests_properties(integrator_stall_threshold PROPERTIES
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} TIMEOUT 120)

add_executable(numprobe_nutrient_pole ${SRC_NOGUI} tests/numerical/nutrient_pole_integrability.c)
target_include_directories(numprobe_nutrient_pole PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(numprobe_nutrient_pole PRIVATE solver_double ${LIBM} ${LIBPTH} ${LIBDL})
add_test(NAME nutrient_pole_integrability COMMAND numprobe_nutrient_pole)
set_tests_properties(nutrient_pole_integrability PROPERTIES
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} TIMEOUT 120)

add_executable(numprobe_state_bound_audit ${SRC_NOGUI} tests/numerical/state_bound_audit.c)
target_include_directories(numprobe_state_bound_audit PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(numprobe_state_bound_audit PRIVATE solver_double ${LIBM} ${LIBPTH} ${LIBDL})
add_test(NAME state_bound_audit COMMAND numprobe_state_bound_audit)
set_tests_properties(state_bound_audit PROPERTIES
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} TIMEOUT 120)


add_executable(single_lsb_sweep ${SRC} tools/single_lsb_sweep.c)
target_include_directories(single_lsb_sweep PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(single_lsb_sweep PRIVATE solver_double PkgConfig::RAYLIB ${LIBM} ${LIBPTH} ${LIBDL} ${LIBGL})

cmake_minimum_required(VERSION 3.20)
project(digital_twin_c C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ---------------- kompilacija ----------------
# TODO: Split compile flags by build type (`Debug` vs `Release`) and avoid forcing `-ffast-math` globally.
 add_compile_options(-Wall -Wextra -O3 -ffast-math)

# ---------------- fajlovi --------------------
add_compile_definitions(SOLVER_REAL_AS=double )
file(GLOB_RECURSE SRC CONFIGURE_DEPENDS src/*.c)
add_executable(simulator ${SRC})
target_include_directories(simulator PRIVATE ${PROJECT_SOURCE_DIR}/include)

# ---------------- solver_double --------------
# 1) nađi statičku ili deljenu biblioteku
# TODO: Replace hardcoded `/usr/local/lib` and `/usr/lib` search paths with a configurable CMake cache variable or package config.
find_library(SOLVER_DOUBLE_PATH
             NAMES solver_double
             PATHS /usr/local/lib /usr/lib
             NO_DEFAULT_PATH
             REQUIRED)

# 2) napravi IMPORTED target
add_library(solver_double STATIC IMPORTED)
set_target_properties(solver_double PROPERTIES
    IMPORTED_LOCATION             ${SOLVER_DOUBLE_PATH}
    INTERFACE_INCLUDE_DIRECTORIES /usr/local/include  # koriguje po potrebi
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

message(STATUS "Solver at : ${SOLVER_DOUBLE_PATH}")
message(STATUS "Sources   : ${SRC}")


# ---------------- Criterion test framework ----------------
pkg_check_modules(CRITERION REQUIRED IMPORTED_TARGET criterion)

file(GLOB_RECURSE TEST_SRC CONFIGURE_DEPENDS tests/*.c)

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


# ---------------- OpenMP -----------------------
# find_package(OpenMP REQUIRED)
#
# if(OpenMP_FOUND)
#     target_link_libraries(simulator PRIVATE OpenMP::OpenMP_C)
#     if(TEST_SRC)
#         target_link_libraries(unit_tests PRIVATE OpenMP::OpenMP_C)
#     endif()
# endif()

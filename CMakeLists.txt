cmake_minimum_required(VERSION 3.20)
project(digital_twin_c C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# ---------------- kompilacija ----------------
# add_compile_options(-Wall -Wextra -O3 -ffast-math)

# ---------------- fajlovi --------------------
add_compile_definitions(SOLVER_REAL_AS=double )
file(GLOB_RECURSE SRC CONFIGURE_DEPENDS src/*.c)
add_executable(simulator ${SRC})
target_include_directories(simulator PRIVATE ${PROJECT_SOURCE_DIR}/include)

# ---------------- solver_double --------------
# 1) nađi statičku ili deljenu biblioteku
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
        solver_double          # ← sada ima punu putanju
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
    list(FILTER SRC EXCLUDE REGEX ".*/main\\.c$")
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
    add_test(NAME all_tests COMMAND unit_tests)
endif()


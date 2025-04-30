cmake_minimum_required(VERSION 3.10)
project(digital_twin_c)

# C standard
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -ffast-math")

# Include folder
include_directories(include)

# Raylib
find_package(raylib REQUIRED)

# Izvori koda
file(GLOB SOURCES "src/*.c" "src/*/*.c")

# Kreiranje izvršne datoteke
add_executable(simulator ${SOURCES})

# Povezivanje sa raylib i njegovim zavisnostima
target_link_libraries(simulator raylib m pthread dl GL)
message(STATUS "Found sources: ${SOURCES}")


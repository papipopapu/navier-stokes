cmake_minimum_required(VERSION 3.10)
project(NavierStokesFluidSimulation VERSION 1.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler warnings
if(MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# Find SFML
find_package(SFML 2.5 COMPONENTS graphics window system REQUIRED)

# Source files
set(SOURCES
    main.cpp
    Grid.cpp
    FluidSolver.cpp
    Renderer.cpp
    Simulation.cpp
)

set(HEADERS
    Config.h
    Grid.h
    FluidSolver.h
    Renderer.h
    Simulation.h
)

# Create executable
add_executable(fluid_simulation ${SOURCES} ${HEADERS})

# Link libraries
target_link_libraries(fluid_simulation 
    PRIVATE 
        sfml-graphics
        sfml-window
        sfml-system
)

# Include directories
target_include_directories(fluid_simulation PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Installation
install(TARGETS fluid_simulation DESTINATION bin)

# Optional: Build examples
option(BUILD_EXAMPLES "Build example simulations" ON)

if(BUILD_EXAMPLES)
    message(STATUS "Building with examples enabled")
endif()

# Print build information
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID}")

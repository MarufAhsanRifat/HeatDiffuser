cmake_minimum_required(VERSION 3.16)
project(HeatDiffuser LANGUAGES Fortran CXX)

# Enable Fortran 2018 features where possible
set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_STANDARD_RQUIRED ON)

# Build the Fortran Library
add_library(heat_fortran STATIC heat_core.f90 heat_grid.f90 heat_grid_oop.f90)

# Build the C++ executablea and link
add_executable(heat_diffuser main.cpp)
target_link_libraries(heat_diffuser PRIVATE heat_fortran)


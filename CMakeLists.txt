cmake_minimum_required(VERSION 3.16)
project(HeatDiffuser LANGUAGES Fortran CXX)

set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

# Default to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()

# Override debug flags (optional)
set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -fbacktrace -fbounds-check" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -fsanitize=address" CACHE STRING "" FORCE)

add_library(heat_fortran STATIC heat_core.f90 heat_grid.f90 heat_grid_oop.f90)
add_executable(heat_diffuser main.cpp)
target_link_libraries(heat_diffuser PRIVATE heat_fortran)

# Optional: add a custom target for quick debug build
add_custom_target(debug
    COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
    COMMAND ${CMAKE_COMMAND} --build .
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Building Debug version"
)

cmake_minimum_required (VERSION 2.8)

project (SPMD_APPLETS)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory (marching_cubes)
add_subdirectory (poisson)
add_subdirectory (contour_tets)


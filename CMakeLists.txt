cmake_minimum_required (VERSION 2.8)

project (SPMD_APPLETS)

set (SCALING_TEST_BUILD OFF CACHE BOOL "build for scaling tests")
if (${SCALING_TEST_BUILD})
  add_definitions("-DSCALING_TEST_BUILD")
endif()

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory (marching_cubes)
add_subdirectory (poisson)
add_subdirectory (contour_tets)


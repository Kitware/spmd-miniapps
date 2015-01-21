cmake_minimum_required (VERSION 2.8)

project (Marching_Cubes)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  "${CMAKE_CURRENT_LIST_DIR}/../cmake/Modules/"
)

find_package (Boost 1.55 REQUIRED chrono system)
find_package (ISPC REQUIRED)

if (CMAKE_SYSTEM_PROCESSOR STREQUAL "k1om")
  set (MIC_BUILD ON)
endif()

set (EXTRA_COMPILER_FLAGS
  "-Wall -Wextra -Wstrict-aliasing -Wno-unused"
  CACHE STRING ""
)
add_definitions (${EXTRA_COMPILER_FLAGS})

set (MARCHING_CUBES_TYPE_DOUBLE ON CACHE BOOL "")
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/type.h.in ${CMAKE_CURRENT_BINARY_DIR}/type.h
)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

set (ISPC_SOURCE_FILES
  "MarchingCubes.ispc"
  "MarchingCubesShortVec.ispc"
  "MarchingCubesV2.ispc"
)

if (MIC_BUILD)
  set (MARCHING_CUBES_ISPC_GANG_SIZE "16" CACHE STRING "")
  ispc_generate_cxx (ISPC_SRCS
    ${MARCHING_CUBES_ISPC_GANG_SIZE}
    "knc-i1x${MARCHING_CUBES_ISPC_GANG_SIZE}.h"
    "-I${CMAKE_CURRENT_BINARY_DIR}"
    ${ISPC_SOURCE_FILES}
  )
  include_directories("${ISPC_ROOT}/examples/intrinsics")
else ()
  set (MARCHING_CUBES_ISPC_GANG_SIZE "16" CACHE STRING "")
  set (MARCHING_CUBES_ISPC_ISA "avx2" CACHE STRING "")
  set (ISPC_TARGET
    "${MARCHING_CUBES_ISPC_ISA}-i32x${MARCHING_CUBES_ISPC_GANG_SIZE}"
  )

  ispc_compile (ISPC_OBJS ISPC_SRCS
    "--target=${ISPC_TARGET} -I${CMAKE_CURRENT_BINARY_DIR}"
    ${ISPC_SOURCE_FILES}
  )
endif (MIC_BUILD)

include_directories (
  ${Boost_INCLUDE_DIRS}
  "${CMAKE_CURRENT_SOURCE_DIR}/../utils"
)

add_executable (MarchingCubes MACOSX_BUNDLE
  MarchingCubes.cxx
  MarchingCubesTables.cxx
  MarchingCubesV2.cxx
  MarchingCubesV2.1.cxx
  main.cxx
  ${ISPC_SRCS}
)

target_link_libraries (MarchingCubes
  ${ISPC_OBJS}
  ${Boost_LIBRARIES}
)


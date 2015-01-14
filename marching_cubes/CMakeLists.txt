cmake_minimum_required (VERSION 2.8)

project (Marching_Cubes)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  "${CMAKE_CURRENT_LIST_DIR}/../cmake/Modules/"
)

find_package (Boost 1.55 REQUIRED chrono system)
find_package (ISPC REQUIRED)

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
  ispc_generate_cxx (ISPC_SRCS
    "16"
    "knc-i1x16.h"
    "-I${CMAKE_CURRENT_BINARY_DIR}"
    ${ISPC_SOURCE_FILES}
  )
  include_directories("${ISPC_ROOT}/examples/intrinsics")
else ()
  ispc_compile (ISPC_OBJS ISPC_SRCS
    "--target=avx2-i32x8 -I${CMAKE_CURRENT_BINARY_DIR}"
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


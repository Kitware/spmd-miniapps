cmake_minimum_required (VERSION 2.8)

project (TetrahedronMeshContour)

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

set (TET_MESH_CONT_TYPE_DOUBLE ON CACHE BOOL "")
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/type.h.in ${CMAKE_CURRENT_BINARY_DIR}/type.h
)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

set (ISPC_SOURCE_FILES
  "Isosurface.ispc"
)

if (MIC_BUILD)
  set (TET_MESH_CONT_ISPC_GANG_SIZE "16" CACHE STRING "")
  set (ISPC_INTRINSICS_HEADER "knc-i1x${TET_MESH_CONT_ISPC_GANG_SIZE}")

  ispc_generate_cxx (ISPC_SRCS
    ${TET_MESH_CONT_ISPC_GANG_SIZE}
    ${ISPC_INTRINSICS_HEADER}
    "-I${CMAKE_CURRENT_BINARY_DIR}"
    ${ISPC_SOURCE_FILES}
  )
  include_directories("${ISPC_ROOT}/examples/intrinsics")
else ()
  set (TET_MESH_CONT_ISPC_GANG_SIZE "16" CACHE STRING "")
  set (TET_MESH_CONT_ISPC_ISA "avx2" CACHE STRING "")
  set (ISPC_TARGET
    "${TET_MESH_CONT_ISPC_ISA}-i32x${TET_MESH_CONT_ISPC_GANG_SIZE}"
  )

  ispc_compile (ISPC_OBJS ISPC_SRCS
    "--target=${ISPC_TARGET} -I${CMAKE_CURRENT_BINARY_DIR}"
    ${ISPC_SOURCE_FILES}
  )
endif (MIC_BUILD)

add_definitions ( -DISPC_GANG_SIZE=${TET_MESH_CONT_ISPC_GANG_SIZE} )

include_directories (
  ${Boost_INCLUDE_DIRS}
  "${CMAKE_CURRENT_SOURCE_DIR}/../utils"
)

add_executable (TetmeshContour MACOSX_BUNDLE
  Isosurface.cxx
  main.cxx
  MarchingTetsTables.cxx
  ${ISPC_SRCS}
)

target_link_libraries (TetmeshContour
  ${ISPC_OBJS}
  ${Boost_LIBRARIES}
)


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
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_COMPILER_FLAGS}")

set (ISPC_SOURCE_FILES
  "Isosurface.ispc"
)

if (MIC_BUILD)
  set (ISPC_GANG_SIZE "16" CACHE STRING "")
  set (ISPC_INTRINSICS_HEADER "knc-i1x${ISPC_GANG_SIZE}")

  ispc_generate_cxx (ISPC_SRCS ${ISPC_GANG_SIZE} ${ISPC_INTRINSICS_HEADER} ""
                     ${ISPC_SOURCE_FILES})
else ()
  set (ISPC_GANG_SIZE "16" CACHE STRING "")
  set (ISPC_ISA "avx2" CACHE STRING "")
  set (ISPC_TARGET "${ISPC_ISA}-i32x${ISPC_GANG_SIZE}")

  ispc_compile (ISPC_OBJS ISPC_SRCS "--target=${ISPC_TARGET}"
                ${ISPC_SOURCE_FILES})
endif (MIC_BUILD)

add_definitions ( -DISPC_GANG_SIZE=${ISPC_GANG_SIZE} )

include_directories (
  ${Boost_INCLUDE_DIRS}
  "${ISPC_ROOT}/examples/intrinsics"
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


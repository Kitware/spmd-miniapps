cmake_minimum_required (VERSION 2.8)

project (Poisson)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  "${CMAKE_CURRENT_LIST_DIR}/../cmake/Modules/"
)

find_package (Boost 1.55 REQUIRED chrono system)
find_package (ISPC REQUIRED)

include_directories (${Boost_INCLUDE_DIRS} SYSTEM)

if (CMAKE_SYSTEM_PROCESSOR STREQUAL "k1om")
  set (MIC_BUILD ON)
endif()

set (EXTRA_COMPILER_FLAGS
  "-Wall -Wextra -Wstrict-aliasing -Wno-unused"
  CACHE STRING ""
)
add_definitions (${EXTRA_COMPILER_FLAGS})

set (POISSON_SOLVER_TYPE_DOUBLE ON CACHE BOOL "")
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/type.h.in ${CMAKE_CURRENT_BINARY_DIR}/type.h
)
include_directories(${CMAKE_CURRENT_BINARY_DIR})


set (ISPC_SOURCE_FILES "poisson.ispc")

if (MIC_BUILD)
  set (POISSON_SOLVER_ISPC_GANG_SIZE "16" CACHE STRING "")
  ispc_generate_cxx (ISPC_SRCS
    ${POISSON_SOLVER_ISPC_GANG_SIZE}
    "knc-i1x${POISSON_SOLVER_ISPC_GANG_SIZE}.h"
    "-I${CMAKE_CURRENT_BINARY_DIR}"
    ${ISPC_SOURCE_FILES}
  )
  include_directories("${ISPC_ROOT}/examples/intrinsics" SYSTEM)
else ()
  set (POISSON_SOLVER_ISPC_GANG_SIZE "16" CACHE STRING "")
  set (POISSON_SOLVER_ISPC_ISA "avx2" CACHE STRING "")
  set (ISPC_TARGET
    "${POISSON_SOLVER_ISPC_ISA}-i32x${POISSON_SOLVER_ISPC_GANG_SIZE}"
  )
  ispc_compile (ISPC_OBJS ISPC_SRCS
    "--target=${ISPC_TARGET} -I${CMAKE_CURRENT_BINARY_DIR}"
    ${ISPC_SOURCE_FILES}
  )
endif (MIC_BUILD)

include_directories ("${CMAKE_CURRENT_SOURCE_DIR}/../utils")

add_executable (poisson MACOSX_BUNDLE
  main.cxx
  poisson.cxx
  ${ISPC_SRCS}
)

target_link_libraries (poisson
  ${ISPC_OBJS}
  ${Boost_LIBRARIES}
)


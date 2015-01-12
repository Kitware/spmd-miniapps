cmake_minimum_required (VERSION 2.8)

project (Marching_Cubes)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  "${CMAKE_CURRENT_LIST_DIR}/../cmake/Modules/"
)

find_package (Boost 1.55 REQUIRED chrono system)
find_package (ISPC REQUIRED)

if (MIC_BUILD)
  ispc_generate_cxx (ISPC_SRCS 4 "knc-i1x16.h" ""
    Gradient.ispc
    MarchingCubesISPC.ispc
  )
else ()
  ispc_compile (ISPC_OBJS ISPC_SRCS "--target=avx2-i32x16"
    poisson.ispc
  )
endif (MIC_BUILD)

include_directories (
  ${Boost_INCLUDE_DIRS}
  "${ISPC_ROOT}/examples/intrinsics"
  "${CMAKE_CURRENT_SOURCE_DIR}/../utils"
)

add_executable (poisson MACOSX_BUNDLE
  main.cxx
  poisson.cxx
  ${ISPC_SRCS}
)

target_link_libraries (poisson
  ${ISPC_OBJS}
  ${Boost_LIBRARIES}
)


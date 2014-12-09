cmake_minimum_required (VERSION 2.8)

project (Marching_Cubes)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  "${CMAKE_CURRENT_LIST_DIR}/../cmake/Modules/"
)

find_package (Boost 1.55 REQUIRED chrono system)
find_package (ISPC REQUIRED)

set (CMAKE_CXX_FLAGS "-msse4.2")

include_directories (
  ${Boost_INCLUDE_DIRS}
  "${ISPC_ROOT}/examples/intrinsics"
)

#ispc_compile (ISPC_OBJS ISPC_SRCS "--target=avx2-i32x8"
#  Gradient.ispc
#  MarchingCubesISPC.ispc
#)

ispc_generate_cxx (ISPC_SRCS 4 "sse4.h" ""
  Gradient.ispc
  MarchingCubesISPC.ispc
)

add_executable (MarchingCubes MACOSX_BUNDLE
  LoadImage3D.cxx
  MarchingCubesISPC.cxx
  MarchingCubesSIMD.cxx
  MarchingCubesSerial.cxx
  MarchingCubesTables.cxx
  main.cxx
  SaveTriangleMesh3D.cxx
  ${ISPC_SRCS}
)

target_link_libraries (MarchingCubes
  ${ISPC_OBJS}
  ${Boost_LIBRARIES}
)


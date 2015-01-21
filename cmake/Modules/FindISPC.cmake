cmake_policy (PUSH)
cmake_minimum_required (VERSION 2.8)
cmake_policy (POP)

message (STATUS "Looking for ISPC...")
find_program (ISPC_COMMAND ispc)
get_filename_component (ISPC_ROOT ${ISPC_COMMAND} PATH CACHE)

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (ISPC DEFAULT_MSG ISPC_COMMAND ISPC_ROOT)

set (ISPC_FLAGS "" CACHE STRING "ISPC compile flags" )
set (ISPC_FLAGS_DEBUG "-g" "-O0" CACHE STRING "ISPC debug compile flags")
set (ISPC_FLAGS_RELEASE "-O1" "-DNDEBUG" CACHE
     STRING "ISPC release compile flags")
set (ISPC_FLAGS_RELWITHDEBINFO "-g" "-O1" "-DNDEBUG" CACHE
     STRING "ISPC release with debug info compile flags")

function (ispc_compile OBJECTS HEADERS OPTIONS)
  separate_arguments (OPTIONS)

  if ("${CMAKE_BUILD_TYPE}" MATCHES "Debug")
    set (ispc_compile_flags ${ISPC_FLAGS} ${ISPC_FLAGS_DEBUG} ${OPTIONS})
  else ("${CMAKE_BUILD_TYPE}" MATCHES "Release")
    set (ispc_compile_flags ${ISPC_FLAGS} ${ISPC_FLAGS_RELEASE} ${OPTIONS})
  else ("${CMAKE_BUILD_TYPE}" MATCHES "RelWithDebInfo")
    set (ispc_compile_flags
         ${ISPC_FLAGS} ${ISPC_FLAGS_RELWITHDEBINFO} ${OPTIONS})
  endif ()

  foreach (filename ${ARGN})
    get_filename_component (base ${filename} NAME)
    set (object ${CMAKE_CURRENT_BINARY_DIR}/${base}.o)
    set (header ${CMAKE_CURRENT_BINARY_DIR}/${base}.h)

    add_custom_command (
      OUTPUT ${object} ${header}
      COMMAND ${ISPC_COMMAND} ${ispc_compile_flags} ${CMAKE_CURRENT_SOURCE_DIR}/${filename} -o ${object} -h ${header}
      DEPENDS ${filename}
      COMMENT "ISPC: Generating object: ${object} and header ${header}"
    )

    set_source_files_properties (${object} GENERATED)
    set_source_files_properties (${header} GENERATED)

    list (APPEND LOBJECTS ${object})
    list (APPEND LHEADERS ${header})
  endforeach ()

  set(${OBJECTS} ${LOBJECTS} PARENT_SCOPE)
  set(${HEADERS} ${LHEADERS} PARENT_SCOPE)
  include_directories (${CMAKE_CURRENT_BINARY_DIR})
endfunction()

function (ispc_generate_cxx SOURCES GANG_SIZE INTRINSICS_HEADER OPTIONS)
  separate_arguments (OPTIONS)

  set (TARGET "generic-x${GANG_SIZE}")

  foreach (filename ${ARGN})
    get_filename_component (base ${filename} NAME)

    set (source ${CMAKE_CURRENT_BINARY_DIR}/${base}.cxx)
    set (header ${CMAKE_CURRENT_BINARY_DIR}/${base}.h)

    set (ispc_compile_flags "--emit-c++" "--target=${TARGET}"
      "--c++-include-file=${INTRINSICS_HEADER}" ${OPTIONS}
    )

    add_custom_command (
      OUTPUT ${source} ${header}
      COMMAND ${ISPC_COMMAND} ${ispc_compile_flags} ${CMAKE_CURRENT_SOURCE_DIR}/${filename} -o ${source} -h ${header}
      DEPENDS ${filename}
      COMMENT "ISPC: Generating cxx: ${object} and header ${header}"
    )

    set_source_files_properties (${source} GENERATED)
    set_source_files_properties (${header} GENERATED)

    list (APPEND LSOURCES ${source} ${header})
  endforeach ()

  set(${SOURCES} ${LSOURCES} PARENT_SCOPE)
  include_directories (${CMAKE_CURRENT_BINARY_DIR})
endfunction()


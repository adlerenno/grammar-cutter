cmake_minimum_required(VERSION 3.10)
project(grammar-cutter VERSION 1.0.0 DESCRIPTION "given a grammar, get an excerpt")

set(CMAKE_C_STANDARD 11)
# Error on undefined references
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
  set(CMAKE_SHARED_LINKER_FLAGS "-Wl")
else()
  set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")
endif()

# Allow support for large files
# add_definitions(-D_LARGEFILE_SOURCE -D_LARGE_FILES -D_FILE_OFFSET_BITS=64)

# Possibility to enable -march=native to increase performance of popcount
option(OPTIMIZE_FOR_NATIVE "Build with -march=native and -mtune=native" ON)
option(TRACE_SYMBOLS "Add the trace symbols for the panic function" ON)
option(CLI "Enable the compilation of the command-line tool" ON)
option(WEB_SERVICE "Enable the compilation of the web-service" OFF)

# configure_file("include/cgraph.cmake.h" "${CMAKE_CURRENT_BINARY_DIR}/cgraph.h" @ONLY)

# get files with "find include src -type f -name '*.c' | sort"
set(SOURCES
  cutter/cutter.h
  cutter/cutter.c
        serp/repair.h
        serp/repair.c
        serp/encoder.h
        serp/encoder.c
        serp/bits.h
        serp/bits.c
        serp/basics.h
)

if(WEB_SERVICE)
  add_definitions(-DWEB_SERVICE)
endif()

set(INCLUDES
  ${CMAKE_CURRENT_BINARY_DIR}
  cutter
  serp
  /usr/local/include
)


add_library(${PROJECT_NAME} SHARED ${SOURCES})
target_link_directories(${PROJECT_NAME} PRIVATE "/usr/local/lib")  # necessary on MAC to find libraries.

target_include_directories(${PROJECT_NAME} PRIVATE ${INCLUDES})

target_link_libraries(${PROJECT_NAME} PRIVATE m) # link with math library

set_target_properties(${PROJECT_NAME} PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER ${CMAKE_CURRENT_BINARY_DIR}/cgraph.h)

# build with -march=native if OPTIMIZE_FOR_NATIVE option set
if(OPTIMIZE_FOR_NATIVE)
  target_compile_options(${PROJECT_NAME} PRIVATE -march=native -mtune=native)
endif()

# adding the trace symbols
if(TRACE_SYMBOLS)
  target_link_options(${PROJECT_NAME} PRIVATE -rdynamic)
endif()

# CLI
if(CLI)
  add_executable(${PROJECT_NAME}-cli main.c)
  add_dependencies(${PROJECT_NAME}-cli ${PROJECT_NAME}) # add library add dependency

  target_include_directories(${PROJECT_NAME}-cli PRIVATE ${INCLUDES})
  target_link_directories(${PROJECT_NAME}-cli PRIVATE "/usr/local/lib")  # necessary on MAC to find libraries.

  target_link_libraries(${PROJECT_NAME}-cli PRIVATE ${PROJECT_NAME})
  if(WEB_SERVICE)
    target_link_libraries(${PROJECT_NAME}-cli PRIVATE microhttpd)
  endif ()


  if(TRACE_SYMBOLS)
    target_link_options(${PROJECT_NAME}-cli PRIVATE -rdynamic)
  endif()
endif()

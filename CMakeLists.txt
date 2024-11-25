cmake_minimum_required(VERSION 3.10)
project(grammar-cutter VERSION 1.0.0 DESCRIPTION "given a grammar, get an excerpt" LANGUAGES C CXX)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
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
        rp/repair.h
        rp/repair.c
        rp/encoder.h
        rp/encoder.c
        rp/bits.h
        rp/bits.c
        rp/basics.h
)

set(INCLUDES
        ${CMAKE_CURRENT_BINARY_DIR}
        cutter
        rp
        /usr/local/include
)

set(SAKAI_SOURCES
        sakai/RePair.hpp
        sakai/RePair.cpp
        sakai/constants/common_functions.hpp
        sakai/constants/constant_numbers.hpp
        sakai/constants/obp_constants.hpp
        sakai/constants/Const.hpp
        sakai/constants/Type.hpp
        sakai/datastructures/bit_vector.cpp
        sakai/datastructures/BitVec.hpp
        sakai/datastructures/BlockVec.hpp
        sakai/datastructures/hash_table_rp.cpp
        sakai/datastructures/hash_table_rp.hpp
        sakai/datastructures/inner_leaf.cpp
        sakai/datastructures/inner_leaf.hpp
        sakai/datastructures/online_bp.cpp
        sakai/datastructures/online_bp.hpp
        sakai/datastructures/RankVec.hpp
        sakai/datastructures/suc_fbtree.cpp
        sakai/datastructures/suc_fbtree.hpp
        sakai/datastructures/succinct_poslp.cpp
        sakai/datastructures/succinct_poslp.hpp
        sakai/datastructures/WBitsVec.hpp
        sakai/datastructures/flc_vector.hpp
        sakai/datastructures/flc_vector.cpp
        sakai/datastructures/out_leaf.hpp
        sakai/datastructures/out_leaf.cpp
        sakai/datastructures/eflc_vector.hpp
        sakai/datastructures/eflc_vector.cpp
        sakai/datastructures/nrs_dic.hpp
        sakai/datastructures/nrs_dic.cpp
        sakai/datastructures/freq_crdict.hpp
        sakai/datastructures/freq_crdict.cpp
        sakai/useful/RSDic.hpp
        sakai/useful/RSDic.cpp
        sakai/useful/Util.hpp
        sakai/useful/permutation.hpp
        sakai/useful/permutation.cpp
        sakai/useful/BitsUtil.hpp
        sakai/useful/BitsUtil.cpp
        sakai/useful/EnumCoder.hpp
        sakai/useful/EnumCoder.cpp
        sakai/useful/MemUtil.hpp
        sakai/useful/ArrayUtil.hpp
        sakai/useful/BasicSearch.hpp
)

set(SAKAI_INCLUDES
        ${CMAKE_CURRENT_SOURCE_DIR}
        sakai
        sakai/constants
        sakai/datastructures
        sakai/useful)

if(WEB_SERVICE)
  add_definitions(-DWEB_SERVICE)
endif()



# Creating C++ Library
add_library(sakai SHARED ${SAKAI_SOURCES})
target_compile_options(sakai PRIVATE -stdlib=libc++ -Ofast -W -Wall -Wno-deprecated -Wc++11-extensions -DNODEBUG)
target_include_directories(sakai PUBLIC ${SAKAI_INCLUDES})

# Creating the C++ connector
add_library(sakai_connector SHARED
        sakai/connector.hpp
        sakai/connector.cpp)
target_compile_options(sakai_connector PRIVATE -stdlib=libc++ -Ofast -W -Wall -Wno-deprecated -Wc++11-extensions -DNODEBUG)
target_include_directories(sakai_connector PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(sakai_connector PRIVATE sakai stdc++)

# Creating main C library
add_library(${PROJECT_NAME} SHARED ${SOURCES})
target_link_directories(${PROJECT_NAME} PRIVATE "/usr/local/lib")  # necessary on MAC to find libraries.

target_include_directories(${PROJECT_NAME} PRIVATE ${INCLUDES})

target_link_libraries(${PROJECT_NAME} PRIVATE m) # link with math library
target_link_libraries(${PROJECT_NAME} PRIVATE sakai_connector stdc++)

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

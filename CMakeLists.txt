
cmake_minimum_required(VERSION 2.8)
PROJECT(MINIPART)

SET(CMAKE_CXX_FLAGS "-g -Wall -std=c++14")
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -Og -Werror -fsanitize=address -fsanitize=undefined -Wno-unknown-pragmas")
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG -DBOOST_UBLAS_NDEBUG -fopenmp")

FIND_PACKAGE(Boost COMPONENTS system filesystem program_options unit_test_framework REQUIRED)

INCLUDE_DIRECTORIES(${MINIPART_SOURCE_DIR}/include)

option(BUILD_TEST "Build unit tests" ON)

SET (SOURCES
  src/hypergraph.cc
  src/io.cc
  src/coarsener.cc
  src/local_search.cc
  src/solver.cc)

add_library(libminipart ${SOURCES})

SET (BIN_SOURCES
  src/main.cc)

ENABLE_TESTING()

IF (BUILD_TEST)
  ADD_SUBDIRECTORY(test)
ENDIF (BUILD_TEST)

add_executable(minipart ${BIN_SOURCES})
target_link_libraries(minipart
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  libminipart
)


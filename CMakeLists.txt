
cmake_minimum_required(VERSION 2.8)
PROJECT(MINIPART)

SET(CMAKE_CXX_FLAGS "-g -Wall -std=c++14")
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -Og -Werror -fsanitize=address -fsanitize=undefined")
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3")

FIND_PACKAGE(Boost COMPONENTS system filesystem program_options unit_test_framework REQUIRED)

INCLUDE_DIRECTORIES(${MINIPART_SOURCE_DIR}/include)

option(BUILD_TEST "Build unit tests" ON)

ENABLE_TESTING()

IF (BUILD_TEST)
  ADD_SUBDIRECTORY(test)
ENDIF (BUILD_TEST)

SET (SOURCES
  src/main.cc)

add_executable(minipart ${SOURCES})
target_link_libraries(minipart
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
)



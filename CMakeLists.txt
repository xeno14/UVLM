project(UVLM)

cmake_minimum_required(VERSION 2.8)

###########################################################
# setting
###########################################################
SET(CMAKE_CXX_FLAGS_RELEASE "-O2")
SET(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
SET(CMAKE_BUILD_TYPE Release)
add_definitions("-Wall -std=c++14")


###########################################################
# dependencies
###########################################################
find_package(gflags REQUIRED)


###########################################################
# testing
###########################################################
enable_testing()
macro(add_gtest test_name src libs)
  set(test_exec ${test_name}.out)
  add_executable(${test_exec} ${src})
  target_link_libraries(${test_exec} gtest gtest_main pthread ${libs})
  add_test(${test_exec} ${CMAKE_CURRENT_BINARY_DIR}/${test_exec})
endmacro(add_gtest)


###########################################################
add_library(morphing src/morphing.cpp)
add_executable(morphing_viewer src/morphing_viewer.cpp)
target_link_libraries(morphing_viewer morphing)

add_library(wing_main src/wing/main.cpp)
target_link_libraries(wing_main gflags pthread)
add_executable(naca00XX src/wing/naca00XX.cpp)
target_link_libraries(naca00XX wing_main)

project(UVLM)

cmake_minimum_required(VERSION 3.0)

###########################################################
# setting
###########################################################
SET(CMAKE_CXX_FLAGS_RELEASE "-O2")
SET(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
SET(CMAKE_BUILD_TYPE Release)
add_definitions("-Wall -std=c++14 -static")
# set(CMAKE_VERBOSE_MAKEFILE 1)

###########################################################
# dependencies
###########################################################
find_package(gflags REQUIRED)
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
find_library(yaml-cpp REQUIRED)



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
# testing
###########################################################

add_definitions("-DUVLM_PROJECT_SOURCE_DIR=\"${PROJECT_SOURCE_DIR}\"")

set(PROTO_DIR "${PROJECT_SOURCE_DIR}/proto")
macro(enable_gflags target)
  target_link_libraries(${target} gflags pthread)
endmacro(enable_gflags)


###########################################################
# targets
###########################################################
add_custom_command(
  OUTPUT ${PROTO_DIR}/uvlm.pb.cc ${PROTO_DIR}/uvlm.pb.h
  WORKING_DIRECTORY ${PROTO_DIR}
  DEPENDS ${PROTO_DIR}/uvlm.proto
  COMMAND "protoc" "uvlm.proto" "--cpp_out=." "--python_out=."
)

add_library(morphing
            src/morphing.cpp
            src/morphing_viewer.cpp
            ${PROTO_DIR}/uvlm.pb.cc)
enable_gflags(morphing)
target_link_libraries(morphing protobuf)
add_executable(morphing_viewer
               src/morphing_viewer_main.cpp)
target_link_libraries(morphing_viewer morphing)

add_executable(result_to_gnuplot
               ${PROTO_DIR}/uvlm.pb.cc
               src/result_to_gnuplot.cpp)
enable_gflags(result_to_gnuplot)
target_link_libraries(result_to_gnuplot protobuf)

add_executable(result_to_vtk
               ${PROTO_DIR}/uvlm.pb.cc
               src/result_to_vtk.cpp)
enable_gflags(result_to_vtk)
target_link_libraries(result_to_vtk protobuf)

add_library(uvlm
            ${PROTO_DIR}/uvlm.pb.cc
            src/shed.cpp
            src/vortex.cpp
            src/uvlm_vortex_ring.cpp
            src/wing_builder.cpp
            src/calc_load/calc_load.cpp
            src/simulator.cpp)
target_link_libraries(uvlm morphing protobuf glog)


add_library(wing
            ${PROTO_DIR}/uvlm.pb.cc
            src/wing/rect.cpp
            src/wing/naca00XX.cpp)
enable_gflags(wing)          
target_link_libraries(wing protobuf)

add_executable(naca00XX src/wing/naca00XX_main.cpp)
enable_gflags(naca00XX)
target_link_libraries(naca00XX wing)


add_subdirectory(test)

add_executable(calc_load src/calc_load/calc_load_main.cpp)
enable_gflags(calc_load)
target_link_libraries(calc_load uvlm glog)

add_executable(flapping_validation src/validation/flapping_validation.cpp)
target_link_libraries(flapping_validation uvlm wing protobuf glog yaml-cpp)

add_executable(plunging_validation src/validation/plunging_validation.cpp)
target_link_libraries(plunging_validation uvlm wing protobuf glog yaml-cpp)

add_executable(fixed_plate src/validation/fixed_plate.cpp)
target_link_libraries(fixed_plate uvlm wing protobuf glog yaml-cpp)

add_executable(2d_velocity src/tools/2d_velocity_main.cpp)
target_link_libraries(2d_velocity uvlm protobuf glog)

add_executable(multiple_flapping_wings src/multiple_flapping_wings_main.cpp)
enable_gflags(multiple_flapping_wings)
target_link_libraries(multiple_flapping_wings uvlm wing protobuf glog yaml-cpp)

###########################################################
# install
###########################################################
install (
  TARGETS
  morphing_viewer
  result_to_gnuplot
  result_to_vtk
  naca00XX
  flapping_validation
  plunging_validation
  fixed_plate
  calc_load
  2d_velocity
  multiple_flapping_wings
  DESTINATION /usr/local/bin)

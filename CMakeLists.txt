project(UVLM)

cmake_minimum_required(VERSION 3.0)

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
set(PROTO_DIR "${PROJECT_SOURCE_DIR}/proto")
message(${PROTO_DIR})

add_custom_command(
  OUTPUT ${PROTO_DIR}/uvlm.pb.cc ${PROTO_DIR}/uvlm.pb.h
  WORKING_DIRECTORY ${PROTO_DIR}
  DEPENDS ${PROTO_DIR}/uvlm.proto
  # COMMAND "protoc" "-I" "${PROTO_DIR}" "${PROTO_DIR}/uvlm.proto" "--cpp_out=." "--python_out=."
  COMMAND "protoc" "uvlm.proto" "--cpp_out=." "--python_out=."
)

add_library(morphing src/morphing.cpp)
add_executable(morphing_viewer
               ${PROTO_DIR}/uvlm.pb.cc
               src/morphing_viewer.cpp)
target_link_libraries(morphing_viewer morphing protobuf gflags pthread)

add_library(uvlm
            uvlm.pb.cc
            src/linear.cpp
            src/shed.cpp
            src/vortex.cpp
            src/uvlm_vortex_ring.cpp)
target_link_libraries(uvlm morphing protobuf)


add_library(wing
            ${PROTO_DIR}/uvlm.pb.cc
            src/wing/naca00XX.cpp)
target_link_libraries(wing protobuf)

add_executable(naca00XX src/wing/naca00XX_main.cpp)
target_link_libraries(naca00XX wing gflags pthread)


add_subdirectory(test)

add_executable(plug src/plug_main.cpp)
target_link_libraries(plug uvlm gflags pthread)

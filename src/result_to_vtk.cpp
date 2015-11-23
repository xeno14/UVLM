/**
 * @file result_to_paraview.cpp
 * @brief Convert snaphot to vtk file
 */


#include "util.h"
#include "vtk/vtk.h"

#include <gflags/gflags.h>
#include <cstdio>
#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>


DEFINE_string(input, "", "input files with format. e.g. \%08d.dat");
DEFINE_string(output, "stdout", "output destination [stdout, path to file]");

using namespace UVLM;

int main(int argc, char* argv[]) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  
  FILE* fp = nullptr;
  if (FLAGS_output == "stdout") {
    fp = stdout;
  } else {
    fp = fopen(FLAGS_output.c_str(), "w");
  }
  CHECK_OPEN(fp);
  std::cerr << "Output: " << FLAGS_output << std::endl;

  std::ifstream ifs(FLAGS_input, std::ios::binary);
  CHECK_OPEN(ifs);
  std::cerr << "Input : " << FLAGS_input << std::endl;

  proto::Snapshot2 snapshot;
  snapshot.ParseFromIstream(&ifs);
  UVLM::vtk::WriteSnapshot2(fp, snapshot);

  return 0;
}

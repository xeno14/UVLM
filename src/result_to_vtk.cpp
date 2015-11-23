/**
 * @file result_to_paraview.cpp
 * @brief Convert snaphot to vtk file
 */


#include "util.h"
#include "vtk/vtk.h"

#include <gflags/gflags.h>
#include <glog/logging.h>
#include <cstdio>
#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>


DEFINE_string(input, "", "input files with format. e.g. \%08d.dat");
DEFINE_string(output, "", "output destination");

using namespace UVLM;

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  google::InstallFailureSignalHandler();
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  FLAGS_logtostderr = true;
  
  FILE* fp = fopen(FLAGS_output.c_str(), "w");
  CHECK(fp);

  std::ifstream ifs(FLAGS_input, std::ios::binary);
  CHECK(ifs);

  LOG(INFO) << FLAGS_input << " --> " << FLAGS_output;
  proto::Snapshot2 snapshot;
  snapshot.ParseFromIstream(&ifs);
  UVLM::vtk::WriteSnapshot2(fp, snapshot);

  return 0;
}

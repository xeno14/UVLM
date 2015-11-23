/**
 * @file result_to_paraview.cpp
 * @brief Convert snaphot to vtk file
 */

#include "util.h"
#include "vtk/vtk.h"

#include <cstdio>
#include <chrono>
#include <fstream>
#include <gflags/gflags.h>
#include <glob.h>
#include <glog/logging.h>
#include <iostream>
#include <vector>

DEFINE_string(input, "", "input files with format. e.g. \%08d.dat");
DEFINE_string(output, "", "output destination");

namespace {

void Writer(const std::string& input, const std::string& output) {
  LOG(INFO) << input << " --> " << output;

  std::ifstream in(input, std::ios::binary);
  CHECK(in);
  UVLM::proto::Snapshot2 snapshot;
  snapshot.ParseFromIstream(&in);
  FILE* out = fopen(output.c_str(), "w");
  CHECK(out);
  UVLM::vtk::WriteSnapshot2(out, snapshot);
}

std::string OutputPath(const std::string& output, const std::size_t index) {
  char path[256];
  sprintf(path, "%s/vtk%08lu.vtk", output.c_str(), index);
  return std::string(path);
}

void GlobWriter(const std::string& input_path, const std::string& output_path) {
  glob_t globbuf;

  int _ = glob(input_path.c_str(), 0, NULL, &globbuf);
  for (uint32_t i = 0; i < globbuf.gl_pathc; i++) {
    Writer(globbuf.gl_pathv[i], OutputPath(output_path, i));
  }
  globfree(&globbuf);
}

}  // anonymous namespace

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  google::InstallFailureSignalHandler();
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  FLAGS_logtostderr = true;

  GlobWriter(FLAGS_input, FLAGS_output);

  return 0;
}

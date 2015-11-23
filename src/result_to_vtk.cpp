/**
 * @file result_to_paraview.cpp
 * @brief Convert snaphot to vtk file
 */

#include "util.h"
#include "vtk/vtk.h"
#include "recordio/recordio.h"

#include <cstdio>
#include <chrono>
#include <fstream>
#include <gflags/gflags.h>
#include <glob.h>
#include <glog/logging.h>
#include <iostream>
#include <vector>

DEFINE_string(input, "", "glob pettern or recordio");
DEFINE_string(output, "", "output destination");
DEFINE_bool(recordio, false, "use recordio instead of glob");
DEFINE_string(prefix, "vtk", "prefix of output vtk files");

namespace {

void Write(const UVLM::proto::Snapshot2& snapshot, const std::string& output) {
  FILE* out;
  CHECK(out = fopen(output.c_str(), "w"));

  UVLM::vtk::WriteSnapshot2(out, snapshot);
}

std::string OutputPath(const std::string& output_path,
                       const std::string& prefix, const std::size_t index) {
  char path[256];
  sprintf(path, "%s/%s%08lu.vtk", output_path.c_str(), prefix.c_str(), index);
  return std::string(path);
}

std::size_t GlobWriter(const std::string& input_path, const std::string& output_path) {
  glob_t globbuf;

  CHECK(glob(input_path.c_str(), 0, NULL, &globbuf) == 0);
  std::size_t index = 0;
  for (index = 0; index < globbuf.gl_pathc; index++) {
    std::ifstream in;
    CHECK((in.open(globbuf.gl_pathv[index], std::ios::binary), in));

    UVLM::proto::Snapshot2 snapshot;
    CHECK(snapshot.ParseFromIstream(&in));

    Write(snapshot, OutputPath(output_path, FLAGS_prefix, index));
  }
  globfree(&globbuf);
  return index;
}

std::size_t RecordioWriter(const std::string& input_path,
                    const std::string& output_path) {
  std::ifstream in;
  CHECK((in.open(input_path, std::ios::binary), in));
  recordio::RecordReader reader(&in);

  std::size_t index = 0;
  UVLM::proto::Snapshot2 snapshot;
  while (reader.ReadProtocolMessage(&snapshot)) {
    Write(snapshot, OutputPath(output_path, FLAGS_prefix, index));
    ++index;
  }
  return index;
}

}  // anonymous namespace

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  google::InstallFailureSignalHandler();
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  FLAGS_logtostderr = true;

  std::size_t num = 0;
  if (FLAGS_recordio) {
    num = RecordioWriter(FLAGS_input, FLAGS_output);
  } else {
    num = GlobWriter(FLAGS_input, FLAGS_output);
  }
  LOG(INFO) << "write " << num << " files with prefix '" << FLAGS_prefix << "'";

  return 0;
}

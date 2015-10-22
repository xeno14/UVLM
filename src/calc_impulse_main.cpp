/**
 * @file calc_impulse_main.cpp
 * @brief Add description here
 */

#include "../proto/uvlm.pb.h"
#include "util.h"
#include "proto_adaptor.h"

#include <fstream>
#include <gflags/gflags.h>
#include <glog/logging.h>

DEFINE_string(input, "", "input Snapshot2");
DEFINE_string(output, "", "output path");


int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InstallFailureSignalHandler();
  FLAGS_logtostderr = true;

  std::ifstream ifs(FLAGS_input, std::ios::binary);
  CHECK((bool)ifs) << "Unable to open " << FLAGS_input;

  UVLM::proto::Snapshot2 snapshot;
  snapshot.ParseFromIstream(&ifs);

  std::vector<UVLM::VortexContainer> containers;
  auto vortices = UVLM::Snapshot2ToContainers(&containers, snapshot);

  std::ofstream ofs(FLAGS_output, std::ios::app);
  const auto wake_offset = CountTotalSize(containers.begin(), containers.end());
  for (auto it = vortices->cbegin() + wake_offset; it != vortices->cend();
       ++it) {
    ofs << it->Impulse() << std::endl;
  }

  return 0;
}

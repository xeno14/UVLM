/**
 * @file dump_circulation.cpp
 * @brief Add description here
 */

#include "../util.h"
#include "../../proto/uvlm.pb.h"
#include "../proto_adaptor.h"

#include <algorithm>
#include <fstream>
#include <gflags/gflags.h>
#include <glog/logging.h>
#include <glob.h>
#include <tuple>

DEFINE_string(input, "", "path to snapshot2");
DEFINE_string(output, "", "output path");


namespace {

void Dump(const std::string& inpath, const std::string& outpath) {
  std::ifstream ifs(inpath, std::ios::binary);
  CHECK((bool)ifs) << "Unable to open " << FLAGS_input;
  std::ofstream ofs(outpath);
  CHECK((bool)ofs) << "Unable to open " << FLAGS_output;

  UVLM::proto::Snapshot2 snapshot;
  snapshot.ParseFromIstream(&ifs);

  std::vector<UVLM::VortexContainer> containers;
  auto vortices = UVLM::Snapshot2ToContainers(&containers, snapshot);

  // -1: wake
  // 0<: wing index
  for (std::size_t index=0; index<containers.size(); ++index) {
    const auto& c = containers[index];
    for (std::size_t i=0; i<c.rows(); i++) {
      for(std::size_t j=0; j<c.cols(); j++) {
        std::vector<double> data;
        data.push_back(c.at(i, j).gamma());
        for (auto node : c.at(i, j).nodes()) {
          data.push_back(node.x());
          data.push_back(node.y());
          data.push_back(node.z());
        }
        std::vector<std::size_t> nums{index, i, j};
        ofs << UVLM::util::join("\t", nums.begin(), nums.end()) << "\t";
        ofs << UVLM::util::join("\t", data.begin(), data.end())
            << std::endl;
      }
      ofs << std::endl;
    }
  }
}

}  // anonymous namespace


int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InstallFailureSignalHandler();
  FLAGS_logtostderr = true;

  Dump(FLAGS_input, FLAGS_output);
  LOG(INFO) << FLAGS_output;

  return 0;
}

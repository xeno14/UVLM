/**
 * @file calc_impulse_main.cpp
 * @brief Add description here
 */

#include "../proto/uvlm.pb.h"
#include "util.h"
#include "proto_adaptor.h"

#include <algorithm>
#include <fstream>
#include <gflags/gflags.h>
#include <glog/logging.h>
#include <glob.h>
#include <tuple>

DEFINE_string(input, "", "path to result dir");
DEFINE_string(output, "", "output path");

std::vector<std::string> GlobResult(const std::string& path) {
  std::vector<std::string> res;
  const std::string pattern = path + "/*";

  glob_t globbuf;
  std::size_t i;

  int ret = glob(pattern.c_str(), 0, NULL, &globbuf);
  CHECK(ret == 0);
  for (i = 0; i < globbuf.gl_pathc; i++) {
    res.emplace_back(globbuf.gl_pathv[i]);
  }
  globfree(&globbuf);
  std::sort(res.begin(), res.end());
  return res;
}

auto CalcImpulse(const std::string& snapshot2_path) {
  std::ifstream ifs(snapshot2_path, std::ios::binary);
  CHECK((bool)ifs) << "Unable to open " << FLAGS_input;

  UVLM::proto::Snapshot2 snapshot;
  snapshot.ParseFromIstream(&ifs);

  const double t = snapshot.t();
  std::vector<UVLM::VortexContainer> containers;
  auto vortices = UVLM::Snapshot2ToContainers(&containers, snapshot);

  Eigen::Vector3d total_impulse;
  const auto wake_offset = CountTotalSize(containers.begin(), containers.end());
  for (auto it = vortices->cbegin() + wake_offset; it != vortices->cend();
       ++it) {
    total_impulse += it->Impulse();
  }
  return std::make_tuple(t, total_impulse);
}

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InstallFailureSignalHandler();
  FLAGS_logtostderr = true;

  auto paths = GlobResult(FLAGS_input);
  std::ofstream ofs(FLAGS_output);
  CHECK((bool)ofs) << "Unable to open " << FLAGS_output;
  for (const auto& path : GlobResult(FLAGS_input)) {
    LOG(INFO) << path;
    Eigen::Vector3d impulse;
    double t;
    std::tie(t, impulse) = CalcImpulse(path);
    std::vector<double> data{t, impulse.x(), impulse.y(), impulse.z()};
    ofs << UVLM::util::join("\t", data.begin(), data.end()) << std::endl;

  }
  return 0;
}

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

namespace {
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

template <class InputIterator>
auto CalcVortexImpulse(InputIterator first, InputIterator last) {
  Eigen::Vector3d res = Eigen::Vector3d::Zero();
  for (auto it = first; it!=last; ++it) {
    res += it->Impulse();
  }
  return res;
}

// TODO think about its deffinition more carafully
template <class InputIterator>
auto CalcPressureImpulse(InputIterator first, InputIterator last) {
  Eigen::Vector3d res;
  for (auto it = first; it!=last; ++it) {
    res += it->Impulse() * (-1);
  }
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

  const auto wakeoffset =
      UVLM::CountTotalSize(containers.begin(), containers.end());
  // TODO multiple wings
  auto vortex_impulse = CalcVortexImpulse(vortices->begin(), vortices->end());
  auto pressure_impulse =
      CalcPressureImpulse(containers[0].begin(), containers[0].end());
  return std::make_tuple(t, vortex_impulse, pressure_impulse);
}
}  // anonymous namespace

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
    Eigen::Vector3d impulse, pressure_impulse;
    double t;
    std::tie(t, impulse, pressure_impulse) = CalcImpulse(path);
    std::vector<double> data{t, impulse.x(), impulse.y(), impulse.z(),
                             pressure_impulse.x(), pressure_impulse.y(),
                             pressure_impulse.z()};
    ofs << UVLM::util::join("\t", data.begin(), data.end()) << std::endl;

  }
  return 0;
}

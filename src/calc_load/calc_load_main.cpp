/**
 * @file main.cpp
 * @brief ２つのSnapshot2から係数を計算する
 */

#include "../../proto/uvlm.pb.h"
#include "../util.h"
#include "../proto_adaptor.h"
#include "calc_load.h"

#include <algorithm>
#include <glob.h>
#include <fstream>
#include <gflags/gflags.h>
#include <glog/logging.h>

DEFINE_string(pattern, "", "pattern for blob");
DEFINE_string(output, "", "output");
DEFINE_double(U, 1.0, "velocity of main stream");
DEFINE_double(alpha, 0, "angle of attack [deg]");

namespace {
using namespace UVLM;

auto ReadSnapshot2(const char* filename) {
  std::ifstream ifs(filename, std::ios::binary);
  CHECK_OPEN(ifs);
  proto::Snapshot2 snapshot2;
  snapshot2.ParseFromIstream(&ifs);
  return snapshot2;
}

auto ReadSnapshot2(const std::vector<std::string>& filenames) {
  std::vector<proto::Snapshot2> res;
  for (const auto& filename : filenames) {
    res.push_back(ReadSnapshot2(filename.c_str()));
  }
  std::sort(res.begin(), res.end(),
            [](const auto& l, const auto& r) { return l.t() < r.t(); });
  return res;
}

auto Snapshot2Paths(const char* pattern) {
    glob_t globbuf;

    int ret = glob(pattern, 0, NULL, &globbuf);
    if (ret == GLOB_NOMATCH) {
      LOG(ERROR) << "Pattern no match";
    }
    
    std::vector<std::string> res;
    for (std::size_t i = 0; i < globbuf.gl_pathc; i++) {
      res.emplace_back(globbuf.gl_pathv[i]);
    }
    globfree(&globbuf);

    return res;
}

/**
 * @param s0 snapshot2 at t-Δt
 * @param s1 snapshot2 at t
 */
auto Calc(const proto::Snapshot2& s0, const proto::Snapshot2& s1) {
  std::vector<UVLM::VortexContainer> c0, c1;
  UVLM::Snapshot2ToContainers(&c0, s0); UVLM::Snapshot2ToContainers(&c1, s1);

  if (c0.size() != c1.size()) {
    LOG(ERROR) << "Container size not match";
  }

  const double dt = s1.t() - s0.t();
  const double RHO = 1.0;
  const double U = FLAGS_U;
  const double ALPHA = FLAGS_alpha * M_PI / 180;
  Eigen::Vector3d Vinfty(U * cos(ALPHA), 0, U * sin(ALPHA));

  std::vector<Eigen::Vector3d> res;

  auto wake_iterator = GetWake(c1);
  for (std::size_t i = 0; i < c1.size(); i++) {
    auto load = CalcLoad(c1[i], c0[i], wake_iterator.first,
                         wake_iterator.second, Vinfty, RHO, dt);
    load /= (0.5 * RHO * U * U * c1[i].chord() * c1[i].span());
    res.emplace_back(load);
  }
  return res;
}

}  // anonymous namespace

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  FLAGS_logtostderr = true;   // TODO 実行時に --logtostderr にすると怒られる

  auto filepaths = Snapshot2Paths(FLAGS_pattern.c_str());
  auto snapshots = ReadSnapshot2(filepaths);

  std::ofstream ofs(FLAGS_output, std::ios::binary);

  for (std::size_t i=1; i<snapshots.size(); i++) {
    LOG(INFO) << "Processing " << i << ": " << filepaths[i - 1] << " & "
              << filepaths[i];
    const double t = snapshots[i].t();
    auto result = Calc(snapshots[i-1], snapshots[i]);
    ofs << t;
    for (const auto& coeff : result) {
      ofs << "\t" << coeff.x() << "\t" << coeff.y() << "\t" << coeff.z();
    }
    ofs << std::endl;
  }

  return 0;
}

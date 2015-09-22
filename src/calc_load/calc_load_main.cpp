/**
 * @file main.cpp
 * @brief ２つのSnapshot2から係数を計算する
 */

#include "../util.h"
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

int main(int argc, char* argv[]) {
  using namespace UVLM::calc_load;
  google::InitGoogleLogging(argv[0]);
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  FLAGS_logtostderr = true;   // TODO 実行時に --logtostderr にすると怒られる

  auto filepaths = Snapshot2Paths(FLAGS_pattern.c_str());
  auto snapshots = UVLM::calc_load::internal::ReadSnapshot2(filepaths);

  std::ofstream ofs(FLAGS_output, std::ios::binary);

  for (std::size_t i=1; i<snapshots.size(); i++) {
    LOG(INFO) << "Processing " << i << ": " << filepaths[i - 1] << " & "
              << filepaths[i];
    const double t = snapshots[i].t();
    auto result = UVLM::calc_load::internal::Calc(snapshots[i-1], snapshots[i]);
    ofs << t;
    for (const auto& coeff : result) {
      ofs << "\t" << coeff.x() << "\t" << coeff.y() << "\t" << coeff.z();
    }
    ofs << std::endl;
  }

  return 0;
}

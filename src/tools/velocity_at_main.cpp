/**
 * @file velocity_at_main.cpp
 * @brief Add description here
 */


#include "../proto_adaptor.h"
#include "../shed.h"
#include "../util.h"
#include "../recordio/recordio_range.h"

#include <cstdio>
#include <fstream>
#include <gflags/gflags.h>
#include <glog/logging.h>
#include <vector>
#include <map>
#include <set>


DEFINE_string(input, "", "snapshot2 recordio");
DEFINE_string(output, "", "output");
DEFINE_double(x, 0, "x");
DEFINE_double(y, 0, "y");
DEFINE_double(z, 0, "z");

namespace {

struct Data {
  double t;
  Eigen::Vector3d v;
};

Data Calc(const UVLM::proto::Snapshot2& snapshot) {
  const double t = snapshot.t();

  std::vector<UVLM::VortexContainer> containers;
  auto vortices = UVLM::Snapshot2ToContainers(&containers, snapshot);

  const Eigen::Vector3d pos{FLAGS_x, FLAGS_y, FLAGS_z};
  Eigen::Vector3d v;
  UVLM::InducedVelocity(&v, pos, vortices->begin(), vortices->end());

  return Data{t, v};
}


}  // anonymous namespace

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InstallFailureSignalHandler();
  FLAGS_logtostderr = true;

  std::ofstream ofs(FLAGS_output);
  CHECK(ofs);

  ofs << "step\tt\tx\ty\tz" << std::endl;

  std::size_t count = 0;
  for (const auto& snapshot :
       recordio::ReaderRange<UVLM::proto::Snapshot2>(FLAGS_input)) {
    auto data = Calc(snapshot);

    ofs << count << "\t"
        << data.t << "\t"
        << data.v.x() << "\t"
        << data.v.y() << "\t"
        << data.v.z() << "\n";
    ++count;
  }
  LOG(INFO) << count << " lines are written."
            << (count == 0 ? " Is input file correct?" : "");

  return 0;
}

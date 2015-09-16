/**
 * @file simulator_main.cpp
 * @brief Add description here
 */

#include "../proto/uvlm.pb.h"
#include "proto_adaptor.h"
#include "simulator.h"
#include "parameter.h"
#include "wing/wing.h"

#include <gflags/gflags.h>
#include <glog/logging.h>

DEFINE_string(input, "", "setting yaml");
DEFINE_string(output, "", "output path");

namespace {
YAML::Node config;
}

auto InitWing() {
  UVLM::proto::Wing wing;

  // 翼の作成
  // NACA0012: AR=8
  DEFINE_PARAM_VERBOSE(int, rows, config["parameter"]);
  DEFINE_PARAM_VERBOSE(int, cols, config["parameter"]);
  UVLM::wing::NACA00XXGenerator(12, 1., 4., PARAM_rows, PARAM_cols)
      .Generate(&wing);
  auto* origin = wing.mutable_origin();
  origin->CopyFrom(UVLM::Vector3dToPoint({0, 0, 0}));
  return wing;
}

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InstallFailureSignalHandler();
  FLAGS_logtostderr = true;

  config = YAML::LoadFile(FLAGS_input);

  UVLM::Morphing m;
  m.set_plug([](double t) { return sin(t); });

  UVLM::proto::Wing wing = InitWing();

  UVLM::simulator::AddWing(wing, m);
  UVLM::simulator::SetInlet(1, 0, 0);
  UVLM::simulator::SetOutputPath(FLAGS_output);

  auto setting = config["setting"];
  DEFINE_PARAM_VERBOSE(int, steps, setting);
  DEFINE_PARAM_VERBOSE(double, dt, setting);
  UVLM::simulator::Start(PARAM_steps, PARAM_dt);

  return 0;
}

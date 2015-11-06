/**
 * @file stanford_baseline_main.cpp
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
DEFINE_string(output_load, "", "output load path");

namespace {
YAML::Node config;
}

auto InitWing(double x0, double y0, double c) {
  const double ar = 6;
  UVLM::proto::Wing wing;
  int rows = 10;
  int cols = 6;
  UVLM::wing::NACA4digitGenerator(83, c, c * ar / 2, rows, cols)
      .Generate(&wing);
  UVLM::wing::SetOrigin(&wing, {x0, y0, 0});
  return wing;
}

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InstallFailureSignalHandler();
  FLAGS_logtostderr = true;

  config = YAML::LoadFile(FLAGS_input);
  DEFINE_PARAM_VERBOSE(double, k, config);
  DEFINE_PARAM_VERBOSE(double, U, config);
  DEFINE_PARAM_VERBOSE(double, chord, config);

  const double K = PARAM_k;
  const double OMEGA = K * 2 * PARAM_U / PARAM_chord;
  const double PHI = M_PI / 4;
  const double alpha = 5.0 / 180.0 * M_PI;
  const double T = 2 * M_PI / OMEGA;
  const double dt = T / 40;
  const std::size_t steps = 50;

  UVLM::Morphing m;
  m.set_flap([&](double t) { return -PHI * cos(OMEGA * (t - 0*dt)); });
  m.set_alpha(alpha);
  UVLM::simulator::AddWing(InitWing(0, 0, PARAM_chord), m);
  UVLM::simulator::SetInlet(PARAM_U, 0, 0);
  UVLM::simulator::SetOutputPath(FLAGS_output);
  UVLM::simulator::SetOutputLoadPath(FLAGS_output_load);

  UVLM::simulator::Start(steps, dt);

  return 0;
}

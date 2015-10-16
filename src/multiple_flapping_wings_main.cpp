/**
 * @file multiple_flapping_wings_main.cpp
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
DEFINE_bool(calc_load, false, "calc load mode");

namespace {
YAML::Node config;
}

auto InitWing(double x0, double y0) {
  UVLM::proto::Wing wing;

  // wing span=1
  // aspect ratio=4
  const auto param = config["parameter"];
  DEFINE_PARAM_VERBOSE(int, rows, param);
  DEFINE_PARAM_VERBOSE(int, cols, param);
  const double span = 1;
  const double chord = span/4;
  UVLM::wing::RectGenerator(chord, span/2, PARAM_rows, PARAM_cols).Generate(&wing);
  UVLM::wing::SetOrigin(&wing, {x0, y0, 0});
  return wing;
}

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InstallFailureSignalHandler();
  FLAGS_logtostderr = true;

  config = YAML::LoadFile(FLAGS_input);
  const auto param = config["parameter"];

  DEFINE_PARAM_VERBOSE(double, U, param);

  const double U = PARAM_U;                           // Upstream velocity
  const double K = 0.5;                               // Reduced frequency
  const double C = 1;                                 // Chord length
  const double OMEGA = 2 * U * K / C;                 // Flapping frequency
  const double PHI = 15 * M_PI / 180;                 // Angle of flapping
  // const double CHORD = wing.chord();
  const double PHI0 = 0;

  UVLM::Morphing m;
  m.set_flap([&](double t) { return PHI * sin(OMEGA * t + PHI0); });
  // TODO position in parameter file
  UVLM::simulator::AddWing(InitWing(0, 0), m);
  UVLM::simulator::AddWing(InitWing(1, 1), m);
  UVLM::simulator::AddWing(InitWing(1, -1), m);
  UVLM::simulator::SetInlet(PARAM_U, 0, 0);
  UVLM::simulator::SetOutputPath(FLAGS_output);
  UVLM::simulator::SetOutputLoadPath(FLAGS_output_load);

  const auto setting = config["setting"];
  DEFINE_PARAM_VERBOSE(int, steps, setting);
  DEFINE_PARAM_VERBOSE(double, dt, setting);

  if (FLAGS_calc_load) {
    // todo
  } else {
    UVLM::simulator::Start(PARAM_steps, PARAM_dt);
  }

  return 0;
}

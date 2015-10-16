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

const double ar = 4;
const double span = 1.2;
const double chord = 1.2 / ar;
const double frequency = 4;  // [Hz]
const double forward_velocity = 15;
}

auto InitWing(double x0, double y0) {
  UVLM::proto::Wing wing;

  // wing span=1
  // aspect ratio=4
  DEFINE_PARAM_VERBOSE(int, rows, config);
  DEFINE_PARAM_VERBOSE(int, cols, config);
  UVLM::wing::RectGenerator(chord, span / 2, PARAM_rows, PARAM_cols)
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

  const double OMEGA = M_PI * 2 * frequency;  // Flapping frequency
  // TODO angle???
  const double PHI = 15 * M_PI / 180;  // Angle of flapping
  const double PHI0 = -M_PI_2;

  UVLM::Morphing m;
  m.set_flap([&](double t) { return PHI * sin(OMEGA * t + PHI0); });
  // TODO position in configeter file
  UVLM::simulator::AddWing(InitWing(0, 0), m);
  UVLM::simulator::AddWing(InitWing(span, span), m);
  UVLM::simulator::AddWing(InitWing(span, -span), m);
  UVLM::simulator::SetInlet(forward_velocity, 0, 0);
  UVLM::simulator::SetOutputPath(FLAGS_output);
  UVLM::simulator::SetOutputLoadPath(FLAGS_output_load);

  DEFINE_PARAM_VERBOSE(int, steps, config);
  DEFINE_PARAM_VERBOSE(double, dt, config);

  if (FLAGS_calc_load) {
    // todo
  } else {
    UVLM::simulator::Start(PARAM_steps, PARAM_dt);
  }

  return 0;
}

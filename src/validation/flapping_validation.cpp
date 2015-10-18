/**
 * @file flapping_validation.cpp
 * @brief Add description here
 */

#include "../../proto/uvlm.pb.h"
#include "../proto_adaptor.h"
#include "../simulator.h"
#include "../parameter.h"
#include "../wing/wing.h"

#include <gflags/gflags.h>
#include <glog/logging.h>

DEFINE_string(input, "", "setting yaml");
DEFINE_string(output, "", "output path");
DEFINE_string(output_load, "", "output load path");
DEFINE_bool(calc_load, false, "calc load mode");
DEFINE_int32(wing_digit, 12, "wing digit for NACA00XX");

namespace {
YAML::Node config;
}

auto InitWing() {
  UVLM::proto::Wing wing;

  // 翼の作成
  // NACA0012: AR=8
  const auto param = config["parameter"];
  DEFINE_PARAM_VERBOSE(int, rows, param);
  DEFINE_PARAM_VERBOSE(int, cols, param);
  UVLM::wing::NACA00XXGenerator(FLAGS_wing_digit, 1., 4., PARAM_rows, PARAM_cols)
      .Generate(&wing);
  UVLM::wing::SetOrigin(&wing, {0, 0, 0});
  return wing;
}

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InstallFailureSignalHandler();
  FLAGS_logtostderr = true;

  config = YAML::LoadFile(FLAGS_input);
  const auto param = config["parameter"];

  UVLM::proto::Wing wing = InitWing();

  DEFINE_PARAM_VERBOSE(double, U, param);
  DEFINE_PARAM_VERBOSE(double, alpha, param);
  DEFINE_PARAM_VERBOSE(double, phi0, param);

  const double U = PARAM_U;                           // Upstream velocity
  const double K = 0.1;                               // Reduced frequency
  const double C = 1;                                 // Chord length
  const double OMEGA = 2 * U * K / C;                 // Flapping frequency
  const double PHI = 15 * M_PI / 180;                 // Angle of flapping
  const double BETA = 4 * M_PI / 180;                 // twising amp at wing tip
  // const double CHORD = wing.chord();
  const double SPAN = wing.span();
  const double ALPHA = PARAM_alpha / 180. * M_PI;
  const double PHI0 = PARAM_phi0;

  UVLM::Morphing m;
  m.set_flap([&](double t) { return PHI * sin(OMEGA * t + PHI0); });
  m.set_twist([&](const Eigen::Vector3d& x0, double t) {
    return BETA * fabs(x0.y()) / SPAN * sin(OMEGA * t + PHI0);
  });
  m.set_alpha(ALPHA);

  UVLM::simulator::AddWing(wing, m);
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

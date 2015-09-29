/**
 * @file simulator_main.cpp
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
  UVLM::wing::NACA00XXGenerator(12, 1., 4., PARAM_rows, PARAM_cols)
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

  const double U = PARAM_U;                           // Upstream velocity
  const double K = 0.1;                               // Reduced frequency
  const double C = 1;                                 // Chord length
  const double OMEGA = 2 * U * K / C;                 // Flapping frequency
  const double PHI = 15 * M_PI / 180;                 // Angle of flapping
  const double BETA = 4 * M_PI / 180;                 // twising amp at wing tip
  // const double CHORD = wing.chord();
  const double SPAN = wing.span();
  const double ALPHA = PARAM_alpha;

  UVLM::Morphing m;
  m.set_flap([OMEGA, PHI](double t) { return PHI * cos(OMEGA * t); });
  m.set_twist([OMEGA, BETA, SPAN](const Eigen::Vector3d& x0, double t) {
    return BETA * fabs(x0.y()) / SPAN * cos(OMEGA * t);
  });

  UVLM::simulator::AddWing(wing, m);
  UVLM::simulator::SetInlet(PARAM_U * cos(ALPHA), 0, PARAM_U * sin(ALPHA));
  UVLM::simulator::SetOutputPath(FLAGS_output);
  UVLM::simulator::SetOutputLoadPath(FLAGS_output_load);

  const auto setting = config["setting"];
  DEFINE_PARAM_VERBOSE(int, steps, setting);
  DEFINE_PARAM_VERBOSE(double, dt, setting);
  UVLM::simulator::Start(PARAM_steps, PARAM_dt);

  return 0;
}

/**
 * @file plunging_validation.cpp
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

namespace {
YAML::Node config;
}

auto InitWing() {
  UVLM::proto::Wing wing;

  // 翼の作成
  // Rect, AR=6
  const auto param = config["parameter"];
  DEFINE_PARAM_VERBOSE(int, rows, param);
  DEFINE_PARAM_VERBOSE(int, cols, param);
  const double chord = 1;
  const double AR = 6;
  const double span = chord * AR / 2;
  UVLM::wing::NACA00XXGenerator(12, chord, span, PARAM_rows, PARAM_cols)
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
  DEFINE_PARAM_VERBOSE(double, k, param);

  const double U = PARAM_U;                           // Upstream velocity
  const double K = PARAM_k;                           // Reduced frequency
  const double CHORD = wing.chord();
  // const double SPAN = wing.span();
  const double OMEGA = 2 * U * K / CHORD;             // Plunging frequency
  const double PLUNG_AMP = 0.175 * CHORD;

  UVLM::Morphing m;
  m.set_plug(
      [OMEGA, PLUNG_AMP](double t) { return PLUNG_AMP * sin(OMEGA * t); });

  UVLM::simulator::AddWing(wing, m);
  UVLM::simulator::SetInlet(U, 0, 0);
  UVLM::simulator::SetOutputPath(FLAGS_output);

  const auto setting = config["setting"];
  DEFINE_PARAM_VERBOSE(int, steps, setting);
  DEFINE_PARAM_VERBOSE(double, dt, setting);
  UVLM::simulator::Start(PARAM_steps, PARAM_dt);

  return 0;
}

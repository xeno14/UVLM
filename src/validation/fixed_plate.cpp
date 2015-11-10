/**
 * @file fixed_plate.cpp
 * @brief Add description here
 */

#include "../uvlm.h"

#include <gflags/gflags.h>
#include <glog/logging.h>

DEFINE_string(input, "", "setting yaml");
DEFINE_string(output, "", "output path");
DEFINE_string(output_load, "", "output load path");
DEFINE_bool(rotate_wing, false, "rotate wing instead of freestream");

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
  DEFINE_PARAM_VERBOSE(double, aspect_ratio, param);
  const double chord = 1;
  const double AR = PARAM_aspect_ratio;
  const double span = chord * AR / 2;
  UVLM::wing::RectGenerator(chord, span, PARAM_rows, PARAM_cols).Generate(&wing);
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

  const double U = PARAM_U;                     // freestream velocity
  const double alpha = DegToRad(PARAM_alpha);   // angle of attack

  UVLM::Morphing m;
  if (FLAGS_rotate_wing) {
    LOG(INFO) << "rotate wing";
    m.set_alpha(alpha);
    UVLM::simulator::SetInlet(U, 0, 0);
  } else {
    LOG(INFO) << "rotate free stream";
    UVLM::simulator::SetInlet(U * cos(alpha), 0, U * sin(alpha));
  }
  UVLM::simulator::AddWing(wing, m);
  UVLM::simulator::SetOutputPath(FLAGS_output);
  UVLM::simulator::SetOutputLoadPath(FLAGS_output_load);

  const auto setting = config["setting"];
  DEFINE_PARAM_VERBOSE(int, steps, setting);
  DEFINE_PARAM_VERBOSE(double, dt, setting);
  UVLM::simulator::Start(PARAM_steps, PARAM_dt);

  return 0;
}

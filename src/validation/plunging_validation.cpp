/**
 * @file plunging_validation.cpp
 * @brief Add description here
 */

#include "../uvlm.h"
#include <fstream>
#include <gflags/gflags.h>
#include <glog/logging.h>

DEFINE_double(k, 0.5, "reduced frequency");
DEFINE_double(h, 0.175, "h/c");
DEFINE_double(alpha, 5, "angle of attack [deg]");
DEFINE_string(result_path, "", "directory to save Snapshot2");
DEFINE_string(load_path, "", "path to aerodynamic loads");
DEFINE_int32(rows, 6, "chordwise num");
DEFINE_int32(cols, 20, "spanwise num");
DEFINE_int32(steps, 50, "number of steps to simulate");
DEFINE_int32(steps_per_cycle, 40, "number of steps per flapping cycle");

using UVLM::simulator::SimpleSimulator;

namespace {
const double AR = 6;
const double CHORD = 1.;
const double SPAN = CHORD * AR;
const double Q = 1;

double omega() {
  return 2 * Q * FLAGS_k / CHORD;
}

void AddWing(SimpleSimulator* simulator) {
  UVLM::Morphing m;
  const double amp = FLAGS_h * CHORD;
  const double o = omega();
  m.set_plug([amp, o](const double t) { return amp * cos(o * t); });
  m.set_alpha(FLAGS_alpha * M_PI / 180.);
  simulator->AddWing(new UVLM::wing::RectGenerator, m, CHORD, SPAN,
                     FLAGS_rows, FLAGS_cols, {0, 0, 0});
}

void Run() {
  SimpleSimulator simulator;
  AddWing(&simulator);
  simulator.set_result_path(FLAGS_result_path);
  simulator.set_load_path(FLAGS_load_path);
  simulator.set_forward_flight({-Q, 0, 0});

  const double dt = 2. * M_PI / omega() / FLAGS_steps_per_cycle;
  simulator.Run(FLAGS_steps, dt);
}
}  // anonymous namespace

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InstallFailureSignalHandler();

  Run();
  return 0;
}

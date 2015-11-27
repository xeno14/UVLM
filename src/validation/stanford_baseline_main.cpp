/**
 * @file stanford_baseline_main.cpp
 * @brief Add description here
 */

#include "../uvlm.h"
#include "../output.h"
#include "../advect_factory.h"
#include <gflags/gflags.h>
#include <glog/logging.h>

DEFINE_string(result_path, "", "directory to save Snapshot2");
DEFINE_string(load_path, "", "path to aerodynamic loads");
DEFINE_int32(rows, 6, "chordwise num");
DEFINE_int32(cols, 20, "spanwise num");
DEFINE_int32(steps, 50, "number of steps to simulate");
DEFINE_int32(steps_per_cycle, 40, "number of steps per flapping cycle");
DEFINE_double(alpha, 5, "angle of attack");
DEFINE_string(advect, "euler", "advection method (see advect_factory.cpp)");

using UVLM::simulator::SimpleSimulator;

namespace {
const double AR = 6;
const double CHORD = 1.;
const double SPAN = CHORD * AR;
const double Kg = 0.1;
const double Q = 1;
const double OMEGA = 2 * Q * Kg / CHORD;

void AddWing(SimpleSimulator* simulator) {
  UVLM::Morphing m;
  const double omega = OMEGA;
  m.set_flap([omega](double t) { return M_PI_4 * cos(omega * t); });
  m.set_alpha(FLAGS_alpha * M_PI / 180.);
  simulator->AddWing(new UVLM::wing::NACA4digitGenerator(83), m, CHORD, SPAN,
                     FLAGS_rows, FLAGS_cols, {0, 0, 0});
}

void Run() {
  SimpleSimulator simulator;
  AddWing(&simulator);
  simulator.set_result_path(FLAGS_result_path);
  simulator.set_load_path(FLAGS_load_path);
  simulator.set_forward_flight({-Q, 0, 0});
  simulator.set_advection(UVLM::advect::AdvectFactory(FLAGS_advect));

  const double dt = 2. * M_PI / OMEGA / FLAGS_steps_per_cycle;
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

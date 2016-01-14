/**
 * @file sudden_acceleration_main.cpp
 * @brief Validation of sudden acceleration of an uncambered rectangular wing
 * into a constant-speed forward flight.
 * For more details, read Katz and Plotkin §13.12 example 1 (p. 429)
 */

#include "../simulator/simple_simulator.h"
#include "../util.h"
#include <gflags/gflags.h>
#include <glog/logging.h>

DEFINE_string(result_path, "", "directory to save Snapshot2");
DEFINE_string(load_path, "", "path to aerodynamic loads");
DEFINE_int32(rows, 4, "chordwise num");
DEFINE_int32(cols, 12, "spanwise num");
DEFINE_int32(steps, 152, "number of steps to simulate");
DEFINE_double(alpha, 5, "angle of attack");
DEFINE_double(dt, 1. / 16., "delta t");
DEFINE_double(AR, 4, "aspect ratio");

using UVLM::simulator::SimpleSimulator;

namespace {

void AddWing(SimpleSimulator* simulator) {
  const double chord = 1;
  const double span = FLAGS_AR * chord;

  UVLM::Morphing m;
  m.set_alpha(UVLM::util::radians(FLAGS_alpha));
  simulator->AddWing(new UVLM::wing::RectGenerator(), m, chord, span,
                     FLAGS_rows, FLAGS_cols, {0, 0, 0});
}

void Run() {
  SimpleSimulator simulator;
  AddWing(&simulator);
  simulator.set_result_path(FLAGS_result_path);
  simulator.set_load_path(FLAGS_load_path);
  simulator.set_forward_flight({-1, 0, 0});

  simulator.Run(FLAGS_steps, FLAGS_dt);
}

}  // anonymous namespace

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InstallFailureSignalHandler();

  Run();
  return 0;
}

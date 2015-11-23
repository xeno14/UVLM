/**
 * @file multiple_flapping_wings_main.cpp
 * @brief Add description here
 */

#include "uvlm.h"
#include "output.h"
#include <fstream>
#include <gflags/gflags.h>
#include <glog/logging.h>

using UVLM::simulator::SimpleSimulator;

DEFINE_string(result_path, "", "directory to save Snapshot2");
DEFINE_string(load_path, "", "path to aerodynamic loads");
DEFINE_int32(rows, 6, "chordwise num");
DEFINE_int32(cols, 20, "spanwise num");
DEFINE_int32(lines, 2, "number of formation lines");
DEFINE_double(xrel, 3, "relative position x");
DEFINE_double(yrel, 6, "relative position y");
DEFINE_double(phase, 0., "phase difference");
DEFINE_int32(steps, 50, "number of steps to simulate");
DEFINE_int32(steps_per_cycle, 40, "number of steps per flapping cycle");
DEFINE_int32(omp_num_threads, 0, "number of threads. if positive, the number is set.");

namespace {
const double AR = 6;
const double ALPHA = 5. * M_PI / 180.;
const double CHORD = 1.;
const double SPAN = CHORD * AR;
const double Kg = 0.1;
const double Q = 1;
const double OMEGA = 2 * Q * Kg / CHORD;

void AddWing(SimpleSimulator* simulator) {
  UVLM::Morphing m_leading;
  const double omega = OMEGA;
  m_leading.set_flap([omega](double t) { return M_PI_4 * cos(omega * t); });
  m_leading.set_alpha(ALPHA);
  simulator->AddWing(m_leading, CHORD, SPAN, FLAGS_rows, FLAGS_cols, {0, 0, 0});
  for (int i = 1; i < FLAGS_lines; i++) {
    UVLM::Morphing m;
    const double dphi = FLAGS_phase * i;
    const double xrel = FLAGS_xrel * i;
    const double yrel = FLAGS_yrel * i;
    m.set_flap(
        [omega, dphi](double t) { return M_PI_4 * cos(omega * t + dphi); });
    m.set_alpha(ALPHA);
    simulator->AddWing(m, CHORD, SPAN, FLAGS_rows, FLAGS_cols, {xrel, yrel, 0});
    simulator->AddWing(m, CHORD, SPAN, FLAGS_rows, FLAGS_cols, {xrel, -yrel, 0});
  }
}

void Run() {
  SimpleSimulator simulator;
  simulator.set_result_path(FLAGS_result_path);
  simulator.set_load_path(FLAGS_load_path);
  simulator.set_forward_flight({-Q, 0, 0});
  AddWing(&simulator);

  const double dt = 2. * M_PI / OMEGA / FLAGS_steps_per_cycle;
  simulator.Run(FLAGS_steps, dt);
}
}  // anonymous namespace

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InstallFailureSignalHandler();

#ifdef _OPENMP
  if (FLAGS_omp_num_threads > 0)
    omp_set_num_threads(FLAGS_omp_num_threads);
#endif

  Run();
  return 0;
}

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

namespace {
double AR = 6;
double ALPHA = 5. * M_PI / 180.;
double CHORD = 1.;
double SPAN = CHORD * AR;
double Kg = 0.1;
double Q = 1;
double OMEGA = 2 * Q * Kg / CHORD;

}  // anonymous namespace

void Run() {
  SimpleSimulator simulator;
  simulator.set_result_path(FLAGS_result_path);

  UVLM::Morphing m;
  const double omega = OMEGA;
  m.set_flap([omega](double t) { return M_PI_4 * cos(omega * t); });
  m.set_alpha(ALPHA);
  simulator.AddWing(m, CHORD, SPAN, 6, 20, {0, 0, 0});
  // simulator.AddWing(m, CHORD, SPAN, 6, 20, {3, 7, 0});
  simulator.set_forward_flight({-Q, 0, 0});

  const double dt =  2. * M_PI / OMEGA / 40;
  simulator.Run(50, dt);
}

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InstallFailureSignalHandler();
  FLAGS_logtostderr = true;

  Run();
  return 0;
}

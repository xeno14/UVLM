/**
 * @file sandwich_main.cpp
 * @brief Add description here
 */

#include "uvlm.h"
#include "output.h"
#include "advect_factory.h"
#include <fstream>
#include <gflags/gflags.h>
#include <glog/logging.h>

using UVLM::simulator::SimpleSimulator;

DEFINE_string(result_path, "", "directory to save Snapshot2");
DEFINE_string(load_path, "", "path to aerodynamic loads");
DEFINE_int32(rows, 6, "chordwise num");
DEFINE_int32(cols, 20, "spanwise num");
DEFINE_int32(lines, 2, "number of formation lines");
DEFINE_double(x1, 3, "relative position of leading to the 2nd.");
DEFINE_double(y1, 6, "relative position of leading to the 2nd.");
DEFINE_double(x2, 6, "relative position of leading to the 3rd.");
DEFINE_double(y2, 12, "relative position of leading to the 3rd.");
DEFINE_double(phase, 0., "phase difference [deg]");
DEFINE_int32(steps, 50, "number of steps to simulate");
DEFINE_int32(steps_per_cycle, 40, "number of steps per flapping cycle");
DEFINE_string(advect, "euler", "advection method (see advect_factory.cpp)");
DEFINE_double(reduced_frequency, 0.1, "reduced frequency");
DEFINE_bool(half, false, "use half of V-shape (i.e. line)");

namespace {
double AR;
int WINGDIGIT;
double ALPHA;
double CHORD;
double SPAN;
double Kg;
double Q;
double OMEGA;

void InitParam() {
  AR = 6;
  WINGDIGIT = 83;
  ALPHA = 5. * M_PI / 180.;
  CHORD = 1.;
  SPAN = CHORD * AR;
  Kg = FLAGS_reduced_frequency;
  Q = 1;
  OMEGA = 2 * Q * Kg / CHORD;
}

void AddWing(SimpleSimulator* simulator) {
  const double omega = OMEGA;

  std::vector<Eigen::Vector3d> origins{
      {0, 0, 0}, {FLAGS_x1, FLAGS_y1, 0}, {FLAGS_x2, FLAGS_y2, 0}};

  double dphi = 0;
  for (const auto& origin : origins) {
    UVLM::Morphing m;
    m.set_flap(
        [omega, dphi](double t) { return -M_PI_4 * cos(omega * t + dphi); });
    m.set_alpha(ALPHA);

    simulator->AddWing(new UVLM::wing::NACA4digitGenerator(WINGDIGIT), m, CHORD,
                       SPAN, FLAGS_rows, FLAGS_cols, origin);

    dphi += FLAGS_phase / 180. * M_PI;
  }
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

  InitParam();
  Run();
  return 0;
}

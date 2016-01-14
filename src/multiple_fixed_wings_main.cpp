/**
 * @file multiple_fixed_wings_main.cpp
 * @brief Add description here
 */

#include "uvlm.h"
#include "output.h"
#include "advect_factory.h"
#include <fstream>
#include <gflags/gflags.h>
#include <glog/logging.h>

using UVLM::simulator::SimpleSimulator;

DEFINE_string(result_path, "", "Snapshot2 recordio");
DEFINE_string(sheet_path, "", "AllVortexSheets recordio");
DEFINE_string(load_path, "", "path to aerodynamic loads");
DEFINE_int32(rows, 6, "chordwise num");
DEFINE_int32(cols, 20, "spanwise num");
DEFINE_int32(lines, 2, "number of formation lines");
DEFINE_double(xrel, 3, "relative position x");
DEFINE_double(yrel, 6, "relative position y");
DEFINE_double(phase, 0., "phase difference [deg]");
DEFINE_int32(steps, 50, "number of steps to simulate");
DEFINE_double(dt, 0.1, "time per step");
DEFINE_string(advect, "euler", "advection method (see advect_factory.cpp)");
DEFINE_bool(half, false, "use half of V-shape (i.e. line)");

namespace {
double AR;
int WINGDIGIT;
double ALPHA;
double CHORD;
double SPAN;
double Q;

void InitParam() {
  AR = 6;
  WINGDIGIT = 83;
  ALPHA = 5. * M_PI / 180.;
  CHORD = 1.;
  SPAN = CHORD * AR;
  Q = 1;
}

void AddWing(SimpleSimulator* simulator) {
  UVLM::Morphing m;
  simulator->AddWing(new UVLM::wing::NACA4digitGenerator(WINGDIGIT), m,
                     CHORD, SPAN, FLAGS_rows, FLAGS_cols, {0, 0, 0});
  for (int i = 1; i < FLAGS_lines; i++) {
    const double xrel = FLAGS_xrel * i;
    const double yrel = FLAGS_yrel * i;
    simulator->AddWing(new UVLM::wing::NACA4digitGenerator(WINGDIGIT), m, CHORD,
                       SPAN, FLAGS_rows, FLAGS_cols, {xrel, yrel, 0});
    if (!FLAGS_half) {
      simulator->AddWing(new UVLM::wing::NACA4digitGenerator(WINGDIGIT), m,
                         CHORD, SPAN, FLAGS_rows, FLAGS_cols, {xrel, -yrel, 0});
    }
  }
}

void Run() {
  SimpleSimulator simulator;
  AddWing(&simulator);
  if (FLAGS_result_path.size()) simulator.set_result_path(FLAGS_result_path);
  if (FLAGS_sheet_path.size()) simulator.set_sheet_path(FLAGS_sheet_path);
  if (FLAGS_load_path.size()) simulator.set_load_path(FLAGS_load_path);
  simulator.set_forward_flight({-Q, 0, 0});
  simulator.set_advection(UVLM::advect::AdvectFactory(FLAGS_advect));

  simulator.Run(FLAGS_steps, FLAGS_dt);
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

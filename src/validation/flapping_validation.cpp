/**
 * @file flapping_validation.cpp
 * @brief Add description here
 * https://books.google.co.jp/books?id=q6oX4o8jSuAC&pg=PA441&hl=ja&source=gbs_toc_r&cad=4#v=onepage&q&f=false
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
DEFINE_double(alpha, 4, "angle of attack");

using UVLM::simulator::SimpleSimulator;

namespace {
const double AR = 8;
const double CHORD = 1.;
const double SPAN = CHORD * AR;
const double Kg = 0.1;
const double Q = 1;
const double OMEGA = 2 * Q * Kg / CHORD;

void AddWing(SimpleSimulator* simulator) {
  const double omega = OMEGA;
  const double span = SPAN;
  const double beta = 4. * M_PI / 180.;
  const double Aphi = 15. * M_PI / 180.;
  simulator->AddWing(
      new UVLM::wing::RectGenerator(),
      UVLM::Morphing()
          .set_alpha(FLAGS_alpha * M_PI / 180.)
          .set_flap([omega, Aphi](double t) { return Aphi * cos(omega * t); })
          .set_twist([omega, beta, span](const Eigen::Vector3d& x0, double t) {
            return beta * fabs(x0.y()) / span * cos(omega * t + M_PI_2);
          }),
      CHORD, SPAN, FLAGS_rows, FLAGS_cols, {0, 0, 0});
}

void Run() {
  SimpleSimulator simulator;
  AddWing(&simulator);
  simulator.set_result_path(FLAGS_result_path);
  simulator.set_load_path(FLAGS_load_path);
  simulator.set_forward_flight({-Q, 0, 0});

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

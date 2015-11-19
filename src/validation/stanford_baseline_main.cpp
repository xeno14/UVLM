/**
 * @file stanford_baseline_main.cpp
 * @brief Add description here
 */

#include "../uvlm.h"
#include "../output.h"
#include "../multiple_sheet/multiple_sheet.h"
#include <fstream>
#include <gflags/gflags.h>
#include <glog/logging.h>

DEFINE_double(Q, 1, "freestream velocity");
DEFINE_double(alpha, 5, "angle of attack [deg]");
DEFINE_double(k, 0.1, "reduced frequency");
DEFINE_double(steps, 50, "number of steps");
DEFINE_string(output, "", "output load file (if empty, use stdout)");
DEFINE_bool(disable_output, false, "disable output of velocities, circulations, etc.");
DEFINE_int32(rows, 6, "rows");
DEFINE_int32(cols, 20, "cols");
DEFINE_string(output_path, "", "directory to save Snapshot2");

using multiple_sheet::MultipleSheet;

namespace {

double AR;
std::size_t ROWS;
std::size_t COLS;
double CHORD;
double SPAN;
double dx;
double dy;
double INF;
double alpha;
double Q;
double Kg;  // reduced frequency
double OMEGA;
double DT;
Eigen::Vector3d U;

std::unique_ptr<std::ostream> load_os;

MultipleSheet<Eigen::Vector3d> wing_pos;
MultipleSheet<Eigen::Vector3d> wing_pos_init;
std::vector<Eigen::Vector3d> wake_pos;
MultipleSheet<double> wing_gamma;
MultipleSheet<double> wing_gamma_prev;
std::vector<double> wake_gamma;
std::vector<Eigen::Vector3d> cpos, cpos_init, normal, tangent;
UVLM::Morphing m;

std::ofstream ofs_morphing("morphing.dat");

void InitParam() {
  AR = 6;
  ROWS = FLAGS_rows;
  COLS = FLAGS_cols;
  CHORD = 1;
  SPAN = CHORD * AR;
  dx = CHORD / ROWS;
  dy = SPAN / COLS;
  INF = 1e20;
  alpha = FLAGS_alpha / 180. * M_PI;
  Q = FLAGS_Q;
  U = Eigen::Vector3d(Q, 0, 0);
  Kg = FLAGS_k;
  OMEGA = 2 * Q * Kg / CHORD;
  wing_gamma.resize(1, ROWS, COLS);
  wing_gamma_prev.resize(1, ROWS, COLS);

  const double o = OMEGA;
  const double phi = 45. / 180. * M_PI;
  m.set_flap([phi, o](double t) { return phi * cos(o * t); });
  m.set_alpha(alpha);
}

std::size_t pos_index(std::size_t i, std::size_t j) {
  return j + i * (COLS + 1);
}
template <class T>
auto& get_pos(T& v, std::size_t i, std::size_t j) {
  return v[j + i * (COLS + 1)];
}
template <class T>
const auto& get_pos(const T& v, std::size_t i, std::size_t j) {
  return v[j + i * (COLS + 1)];
}

std::size_t panel_index(std::size_t i, std::size_t j) { return j + i * COLS; }
template <class T>
auto& get_panel(T& v, std::size_t i, std::size_t j) {
  // CHECK(v.size() == (COLS)*(ROWS));
  return v[j + i * COLS];
}
template <class T>
const auto& get_panel(const T& v, std::size_t i, std::size_t j) {
  // CHECK(v.size() == (COLS)*(ROWS));
  return v[j + i * COLS];
}

void InitPosition(MultipleSheet<Eigen::Vector3d>& pos) {
  UVLM::proto::Wing wing, half;
  UVLM::wing::NACA4digitGenerator(83, CHORD, SPAN / 2, ROWS, COLS / 2)
      .Generate(&half);
  UVLM::wing::WholeWing(&wing, half);
  auto points = UVLM::PointsToVector(wing.points());
  pos.resize(1, ROWS + 1, COLS + 1);
  std::copy(points.begin(), points.end(), pos.begin());
}

auto CollocationPoints(const MultipleSheet<Eigen::Vector3d>& pos) {
  std::vector<Eigen::Vector3d> res;
  for (std::size_t i = 0; i < pos.rows() - 1; i++) {
    for (std::size_t j = 0; j < pos.cols() - 1; j++) {
      res.push_back(
          (pos.at(0, i, j) + pos.at(0, i+1, j) + pos.at(0, i+1, j+1) +
           pos.at(0, i, j+1)) / 4);
    }
  }
  return res;
}

auto Normals(const MultipleSheet<Eigen::Vector3d>& pos) {
  std::vector<Eigen::Vector3d> res;
  for (std::size_t i = 0; i < ROWS; ++i) {
    for (std::size_t j = 0; j < COLS; ++j) {
      Eigen::Vector3d n = (pos.at(0, i + 1, j + 1) - pos.at(0, i, j))
                              .cross(pos.at(0, i, j + 1) - pos.at(0, i + 1, j));
      n.normalize();
      res.emplace_back(n);
    }
  }
  return res;
}

auto VORTEX(const Eigen::Vector3d& x, const Eigen::Vector3d& x1,
            const Eigen::Vector3d& x2, double gamma) {
  // (10.16)
  // impl p. 584
  Eigen::Vector3d res;
  UVLM::BiotSavartLaw(&res, x1, x2, x);
  res *= gamma;
  return res;
}

template <class InputIterator>
Eigen::Vector3d VORING(const Eigen::Vector3d& x, 
    const MultipleSheet<Eigen::Vector3d>& pos,
                       InputIterator gamma_first, std::size_t i,
                       std::size_t j) {
  Eigen::Vector3d u = Eigen::Vector3d::Zero();
  double gamma = *(gamma_first + panel_index(i, j));
  auto p0 = pos.at(0, i, j);
  auto p1 = pos.at(0, i, j + 1);
  auto p2 = pos.at(0, i + 1, j + 1);
  auto p3 = pos.at(0, i + 1, j);
  u += VORTEX(x, p0, p1, gamma);
  u += VORTEX(x, p1, p2, gamma);
  u += VORTEX(x, p2, p3, gamma);
  u += VORTEX(x, p3, p0, gamma);
  return u;
}

template <class InputIterator1, class InputIterator2>
Eigen::Vector3d VORING(const Eigen::Vector3d& x, InputIterator1 pos_first,
                       InputIterator2 gamma_first, std::size_t i,
                       std::size_t j) {
  Eigen::Vector3d u = Eigen::Vector3d::Zero();
  double gamma = *(gamma_first + panel_index(i, j));
  auto p0 = pos_first + pos_index(i, j);
  auto p1 = pos_first + pos_index(i, j + 1);
  auto p2 = pos_first + pos_index(i + 1, j + 1);
  auto p3 = pos_first + pos_index(i + 1, j);
  u += VORTEX(x, *p0, *p1, gamma);
  u += VORTEX(x, *p1, *p2, gamma);
  u += VORTEX(x, *p2, *p3, gamma);
  u += VORTEX(x, *p3, *p0, gamma);
  return u;
}

Eigen::Vector3d BoundVelocity(const Eigen::Vector3d& x) {
  Eigen::Vector3d res = Eigen::Vector3d::Zero();
  for (std::size_t i = 0; i < ROWS; ++i) {
    for (std::size_t j = 0; j < COLS; ++j) {
      res += VORING(x, wing_pos, wing_gamma.begin(), i, j);
    }
  }
  return res;
}

Eigen::Vector3d WakeVelocity(const Eigen::Vector3d& x) {
  Eigen::Vector3d res = Eigen::Vector3d::Zero();
  std::size_t rows = wake_pos.size() / (COLS + 1);
  for (std::size_t i = 0; i < rows - 1; ++i) {
    for (std::size_t j = 0; j < COLS; ++j) {
      res += VORING(x, wake_pos.cbegin(), wake_gamma.cbegin(), i, j);
    }
  }
  return res;
}

auto CalcMatrix() {
  // A_kl
  Eigen::MatrixXd res(ROWS * COLS, ROWS * COLS);
  std::vector<double> gamma(wing_gamma.size(), 1);

  // loop for all bound vortices
  for (std::size_t i = 0; i < ROWS; ++i) {
    for (std::size_t j = 0; j < COLS; ++j) {
      const auto k = panel_index(i, j);
      const auto cp = get_panel(cpos, i, j);
      const auto n = get_panel(normal, i, j);

      // loop for vortex ring
      for (std::size_t ii = 0; ii < ROWS; ++ii) {
        for (std::size_t jj = 0; jj < COLS; ++jj) {
          const auto l = panel_index(ii, jj);
          auto u = VORING(cp, wing_pos.begin(), gamma.cbegin(), ii, jj);
          res(k, l) = u.dot(n);
        }
      }
    }
  }
  return res;
}

Eigen::Vector3d Velocity(const Eigen::Vector3d& x) {
  return U + BoundVelocity(x) + WakeVelocity(x);
}

Eigen::Vector3d MorphingVelocity(const Eigen::Vector3d& x0, double t) {
  Eigen::Vector3d res = Eigen::Vector3d::Zero();
  m.Velocity(&res, x0, t);
  return res;
}

auto CalcRhs(double t) {
  const std::size_t sz = cpos.size();
  Eigen::VectorXd res(sz);
  for (std::size_t K = 0; K < sz; ++K) {
    Eigen::Vector3d u =
        U + WakeVelocity(cpos[K]) - MorphingVelocity(cpos_init[K], t);
    res(K) = -u.dot(normal[K]);
  }
  return res;
}

struct VortexLine {
  Eigen::Vector3d p0, p1;
  Eigen::Vector3d p0_init, p1_init;
  double g;
};

std::vector<VortexLine> GetLines() {
  std::vector<VortexLine> res;
  for (std::size_t i = 0; i < ROWS; i++) {
    for (std::size_t j = 0; j < COLS; j++) {
      std::vector<Eigen::Vector3d> corner = {
          wing_pos.at(0, i, j), wing_pos.at(0, i, j + 1),
          wing_pos.at(0, i + 1, j + 1), wing_pos.at(0, i + 1, j)};
      std::vector<Eigen::Vector3d> corner_init = {
          wing_pos_init.at(0, i, j), wing_pos_init.at(0, i, j + 1),
          wing_pos_init.at(0, i + 1, j + 1), wing_pos_init.at(0, i + 1, j)};
      for (std::size_t k = 0; k < corner.size(); k++) {
        if (i == ROWS - 1 && k == 2) continue;  // skip T.E
        res.push_back(VortexLine{corner[k], corner[(k + 1) % corner.size()],
                                 corner_init[k],
                                 corner_init[(k + 1) % corner_init.size()],
                                 wing_gamma.at(0, i, j)});
      }
    }
  }
  return res;
}

// Simpson's method
Eigen::Vector3d CalcLift2(double t) {
  // steady part
  auto lines = GetLines();
  double Fx = 0, Fy = 0, Fz = 0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : Fx, Fy, Fz)
#endif
  for (std::size_t i = 0; i < lines.size(); i++) {
    const auto& line = lines[i];
    Eigen::Vector3d mp = (line.p0 + line.p1) / 2;
    Eigen::Vector3d mp_init = (line.p0_init + line.p1_init) / 2;
    Eigen::Vector3d u = Velocity(mp) - MorphingVelocity(mp_init, t);
    Eigen::Vector3d df = u.cross(line.p1 - line.p0) * line.g;
    Fx += df.x();
    Fy += df.y();
    Fz += df.z();
  }
  return Eigen::Vector3d(Fx, Fy, Fz);
}

Eigen::Vector3d CalcLift2_unst() {
  // unsteady part
  double Fx = 0, Fy = 0, Fz = 0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : Fx, Fy, Fz)
#endif
  for (std::size_t i = 0; i < ROWS; i++) {
    for (std::size_t j = 0; j < COLS; j++) {
      const double A =
          ((wing_pos.at(0, i+1, j) - wing_pos.at(0, i,j)).cross(
            wing_pos.at(0, i,j+1) - wing_pos.at(0, i, j))).norm();
      const double dg_dt =
          (wing_gamma.at(0, i, j) - wing_gamma_prev.at(0, i, j)) / DT;
      Eigen::Vector3d df = get_panel(normal, i, j) * dg_dt * A;
      Fx += df.x();
      Fy += df.y();
      Fz += df.z();
    }
  }
  return Eigen::Vector3d(Fx, Fy, Fz);
}

void Output(std::size_t step) {
  char filename[256];
  UVLM::proto::Snapshot2 snapshot;
  UVLM::output::SimpleAppendSnapshot(&snapshot, wing_pos.begin(),
                                     wing_pos.end(), wing_gamma.begin(),
                                     wing_gamma.end(), COLS);
  if (wake_gamma.size()) {
    UVLM::output::SimpleAppendSnapshot(&snapshot, wake_pos.cbegin(),
                                       wake_pos.cend(), wake_gamma.cbegin(),
                                       wake_gamma.cend(), COLS);
  }

  sprintf(filename, "%s/%08lu", FLAGS_output_path.c_str(), step);
  std::ofstream ofs(filename);
  CHECK(ofs);
  snapshot.SerializeToOstream(&ofs);

  for (std::size_t K=0; K<wing_pos.size(); ++K) {
    auto u = MorphingVelocity(wing_pos_init[K], step * DT);
    ofs_morphing << wing_pos[K].transpose() << " " << u.transpose() << std::endl;
  }
  ofs_morphing << std::endl << std::endl;
}

void MainLoop(std::size_t step) {
  const double t = step * DT;

  std::copy(wing_gamma.begin(), wing_gamma.end(), wing_gamma_prev.begin());

  // TODO morphing
  CHECK(wing_pos.size() == wing_pos_init.size());
  for (std::size_t i = 0; i < wing_pos.size(); i++) {
    m.Perfome(&wing_pos[i], wing_pos_init[i], t);
  }

  // shed wake
  LOG(INFO) << "Shed";
  std::vector<Eigen::Vector3d> new_wake_pos;
  std::vector<double> new_wake_gamma;
  for (std::size_t j = 0; j <= COLS; j++) {
    new_wake_pos.push_back(wing_pos.at(0 , ROWS, j));
  }
  if (step > 1) {
    for (std::size_t j = 0; j < COLS; j++) {
      new_wake_gamma.push_back(wing_gamma.at(0, ROWS - 1, j));
    }
  }
  new_wake_pos.insert(new_wake_pos.end(), wake_pos.begin(), wake_pos.end());
  new_wake_gamma.insert(new_wake_gamma.end(), wake_gamma.begin(),
                        wake_gamma.end());
  wake_pos.swap(new_wake_pos);
  wake_gamma.swap(new_wake_gamma);

  cpos = CollocationPoints(wing_pos);
  normal = Normals(wing_pos);

  // solve linear
  LOG(INFO) << "Linear";
  auto A = CalcMatrix();
  auto rhs = CalcRhs(t);
  Eigen::FullPivLU<Eigen::MatrixXd> solver(A);
  Eigen::VectorXd gamma_v = solver.solve(rhs);
  for (std::size_t K = 0; K < ROWS * COLS; ++K) wing_gamma[K] = gamma_v(K);

  // calc load
  LOG(INFO) << "Load: joukowski";
  const auto F = CalcLift2(t) + CalcLift2_unst();
  LOG(INFO) << F.transpose();
  const auto C = F / (0.5 * Q * Q * CHORD * SPAN);
  if (FLAGS_output.size()) {
    *load_os << step* DT* Q / CHORD << " " << C.x() << " " << C.z()
             << std::endl;
  } else {
    std::cout << step* DT* Q / CHORD << " " << C.x() << " " << C.z()
              << std::endl;
  }

  // output
  if (!FLAGS_disable_output) Output(step);

  if (step==FLAGS_steps) return;

  // advection
  LOG(INFO) << "Advect";
  std::vector<Eigen::Vector3d> wake_vel(wake_pos.size());

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (std::size_t i = 0; i < wake_pos.size(); i++) {
    wake_vel[i] = Velocity(wake_pos[i]);
  }

  for (std::size_t j = 0; j <= COLS; j++) {
    wake_vel[j] = U;
  }
  for (std::size_t i = 0; i < wake_pos.size(); i++) {
    wake_pos[i] += wake_vel[i] * DT;
  }
}

void SimulatorBody() {
  InitPosition(wing_pos_init);
  wing_pos.resize(wing_pos_init.num(), wing_pos_init.rows(),
                  wing_pos_init.cols());
  std::copy(wing_pos_init.begin(), wing_pos_init.end(), wing_pos.begin());
  cpos_init = CollocationPoints(wing_pos_init);
  DT = 2 * M_PI / OMEGA / 40;
  for (std::size_t i = 1; i <= FLAGS_steps; i++) {
    LOG(INFO) << "step=" << i;
    MainLoop(i);
  }
}

}  // anonymous namespace

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InstallFailureSignalHandler();
  FLAGS_logtostderr = true;
  load_os = std::unique_ptr<std::ostream>(new std::ofstream(FLAGS_output));
  InitParam();
  SimulatorBody();
  return 0;
}

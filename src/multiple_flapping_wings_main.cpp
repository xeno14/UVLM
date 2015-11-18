/**
 * @file multiple_flapping_wings_main.cpp
 * @brief Add description here
 */

#include "uvlm.h"
#include "output.h"
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

namespace {

double AR;
std::size_t ROWS;
std::size_t COLS;
std::size_t NUM;    // number of wings
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

std::vector<Eigen::Vector3d> wing_pos;
std::vector<Eigen::Vector3d> wing_pos_init;
std::vector<Eigen::Vector3d> wake_pos;
std::vector<double> wing_gamma;
std::vector<double> wing_gamma_prev;
std::vector<double> wake_gamma;
std::vector<Eigen::Vector3d> cpos, cpos_init, normal, tangent;
UVLM::Morphing m;

std::ofstream ofs_morphing("morphing.dat");

void InitParam() {
  AR = 6;
  NUM = 2;
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
  wing_gamma.resize(NUM * ROWS * COLS, 0);

  const double o = OMEGA;
  const double phi = 45. / 180. * M_PI;
  m.set_flap([phi, o](double t) { return phi * cos(o * t); });
  m.set_alpha(alpha);
}

template <class T>
auto& get_pos(T& v, std::size_t i, std::size_t j) {
  return v[j + i * (COLS + 1)];
}
template <class T>
const auto& get_pos(const T& v, std::size_t i, std::size_t j) {
  return v[j + i * (COLS + 1)];
}
template <class T>
auto& get_pos(T& v, std::size_t n, std::size_t i, std::size_t j) {
  return v[j + i * (COLS + 1) + n * (ROWS + 1) * (COLS + 1)];
}
template <class T>
const auto& get_pos(const T& v, std::size_t n, std::size_t i, std::size_t j) {
  return v[j + i * (COLS + 1) + n * (ROWS + 1) * (COLS + 1)];
}

std::size_t panel_index(std::size_t i, std::size_t j) { return j + i * COLS; }
template <class T>
auto& get_panel(T& v, std::size_t i, std::size_t j) {
  return v[j + i * COLS];
}
template <class T>
const auto& get_panel(const T& v, std::size_t i, std::size_t j) {
  return v[j + i * COLS];
}
template <class T>
auto& get_panel(T& v, std::size_t n, std::size_t i, std::size_t j) {
  return v[j + i * COLS + n * ROWS * COLS];
}
template <class T>
const auto& get_panel(const T& v, std::size_t n, std::size_t i, std::size_t j) {
  return v[j + i * COLS + n * ROWS * COLS];
}

template <class T>
auto pos_cbegin(const T& v, std::size_t n) {
  return v.cbegin() + n * (ROWS + 1) * (COLS + 1);
}
template <class T>
auto pos_cend(const T& v, std::size_t n) {
  return v.cbegin() + (n + 1) * (ROWS + 1) * (COLS + 1);
}
template <class T>
auto panel_cbegin(const T& v, std::size_t n) {
  return v.cbegin() + n * ROWS * COLS;
}
template <class T>
auto panel_cend(const T& v, std::size_t n) {
  return v.cbegin() + (n + 1) * ROWS * COLS;
}

auto InitPosition(const std::vector<Eigen::Vector3d>& origins) {
  CHECK(NUM == origins.size()) << "size mismatch";
  std::vector<Eigen::Vector3d> res;
  for (auto origin : origins) {
    UVLM::proto::Wing wing, half;
    UVLM::wing::NACA4digitGenerator(83, CHORD, SPAN / 2, ROWS, COLS / 2)
        .Generate(&half);
    UVLM::wing::WholeWing(&wing, half);
    auto pos = UVLM::PointsToVector(wing.points());
    std::for_each(pos.begin(), pos.end(), [origin](auto& p) { p += origin; });
    res.insert(res.end(), pos.begin(), pos.end());
  }
  return res;
}

auto CollocationPoints(const std::vector<Eigen::Vector3d>& pos) {
  std::vector<Eigen::Vector3d> res(ROWS * COLS);
  for (std::size_t i = 0; i < ROWS; i++) {
    for (std::size_t j = 0; j < COLS; j++) {
      auto& p = get_panel(res, i, j);
      p = (get_pos(pos, i, j) + get_pos(pos, i + 1, j) +
           get_pos(pos, i + 1, j + 1) + get_pos(pos, i, j + 1)) /
          4;
    }
  }
  return res;
}

auto Normals(const std::vector<Eigen::Vector3d>& pos) {
  std::vector<Eigen::Vector3d> res(ROWS * COLS);
  for (std::size_t i = 0; i < ROWS; ++i) {
    for (std::size_t j = 0; j < COLS; ++j) {
      Eigen::Vector3d n =
          (get_pos(pos, i + 1, j + 1) - get_pos(pos, i, j))
              .cross(get_pos(pos, i, j + 1) - get_pos(pos, i + 1, j));
      n.normalize();
      get_panel(res, i, j) = n;
    }
  }
  return res;
}

auto Tangents(const std::vector<Eigen::Vector3d>& pos) {
  std::vector<Eigen::Vector3d> res(ROWS * COLS);
  for (std::size_t i = 0; i < ROWS; ++i) {
    for (std::size_t j = 0; j < COLS; ++j) {
      Eigen::Vector3d t = ((get_pos(pos, i + 1, j + 1) - get_pos(pos, i, j)) +
                           (get_pos(pos, i + 1, j) - get_pos(pos, i, j + 1))) /
                          2;
      t.normalize();
      get_panel(res, i, j) = t;
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

Eigen::Vector3d VORING(const Eigen::Vector3d& x,
                       const std::vector<Eigen::Vector3d>& pos,
                       const std::vector<double>& gammas, std::size_t i,
                       std::size_t j) {
  Eigen::Vector3d u = Eigen::Vector3d::Zero();
  double gamma = get_panel(gammas, i, j);
  u += VORTEX(x, get_pos(pos, i, j), get_pos(pos, i, j + 1), gamma);
  u += VORTEX(x, get_pos(pos, i, j + 1), get_pos(pos, i + 1, j + 1), gamma);
  u += VORTEX(x, get_pos(pos, i + 1, j + 1), get_pos(pos, i + 1, j), gamma);
  u += VORTEX(x, get_pos(pos, i + 1, j), get_pos(pos, i, j), gamma);
  return u;
}

Eigen::Vector3d BoundVelocity(const Eigen::Vector3d& x) {
  Eigen::Vector3d res = Eigen::Vector3d::Zero();
  for (std::size_t i = 0; i < ROWS; ++i) {
    for (std::size_t j = 0; j < COLS; ++j) {
      res += VORING(x, wing_pos, wing_gamma, i, j);
    }
  }
  return res;
}

Eigen::Vector3d WakeVelocity(const Eigen::Vector3d& x) {
  Eigen::Vector3d res = Eigen::Vector3d::Zero();
  std::size_t rows = wake_pos.size() / (COLS + 1);
  for (std::size_t i = 0; i < rows - 1; ++i) {
    for (std::size_t j = 0; j < COLS; ++j) {
      res += VORING(x, wake_pos, wake_gamma, i, j);
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
          auto u = VORING(cp, wing_pos, gamma, ii, jj);
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
          get_pos(wing_pos, i, j), get_pos(wing_pos, i, j + 1),
          get_pos(wing_pos, i + 1, j + 1), get_pos(wing_pos, i + 1, j)};
      std::vector<Eigen::Vector3d> corner_init = {
          get_pos(wing_pos_init, i, j), get_pos(wing_pos_init, i, j + 1),
          get_pos(wing_pos_init, i + 1, j + 1),
          get_pos(wing_pos_init, i + 1, j)};
      for (std::size_t k = 0; k < corner.size(); k++) {
        if (i == ROWS - 1 && k == 2) continue;  // skip T.E
        res.push_back(VortexLine{corner[k], corner[(k + 1) % corner.size()],
                                 corner_init[k],
                                 corner_init[(k + 1) % corner_init.size()],
                                 get_panel(wing_gamma, i, j)});
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
          ((get_pos(wing_pos, i + 1, j) - get_pos(wing_pos, i, j))
               .cross(get_pos(wing_pos, i, j + 1) - get_pos(wing_pos, i, j)))
              .norm();
      const double dg_dt =
          (get_panel(wing_gamma, i, j) - get_panel(wing_gamma_prev, i, j)) / DT;
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
  for (std::size_t n=0; n<NUM; n++) {
    UVLM::output::SimpleAppendSnapshot(&snapshot,
        pos_cbegin(wing_pos, n), pos_cend(wing_pos, n),
        panel_cbegin(wing_gamma, n), panel_cend(wing_gamma, n),
        COLS);
  }
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

  wing_gamma_prev = wing_gamma;

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
    new_wake_pos.push_back(get_pos(wing_pos, ROWS, j));
  }
  if (step > 1) {
    for (std::size_t j = 0; j < COLS; j++) {
      new_wake_gamma.push_back(get_panel(wing_gamma, ROWS - 1, j));
    }
  }
  new_wake_pos.insert(new_wake_pos.end(), wake_pos.begin(), wake_pos.end());
  new_wake_gamma.insert(new_wake_gamma.end(), wake_gamma.begin(),
                        wake_gamma.end());
  wake_pos.swap(new_wake_pos);
  wake_gamma.swap(new_wake_gamma);

  cpos = CollocationPoints(wing_pos);
  normal = Normals(wing_pos);
  tangent = Tangents(wing_pos);

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
  wing_pos_init = InitPosition({{0, 0, 0}, {2*CHORD, 1.5*SPAN, 0}});
  wing_pos = wing_pos_init;
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

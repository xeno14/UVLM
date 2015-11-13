/**
 * @file unsteady_plate.cpp
 * @brief Add description here
 */

#include "../uvlm.h"
#include <fstream>
#include <gflags/gflags.h>
#include <glog/logging.h>

DEFINE_double(AR, 4, "aspect ratio");
DEFINE_double(Q, 1, "freestream velocity");
DEFINE_double(alpha, 5, "angle of attack [deg]");

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
Eigen::Vector3d U;


std::vector<Eigen::Vector3d> wing_pos;
std::vector<Eigen::Vector3d> wake_pos;
std::vector<double> wing_gamma;
std::vector<double> wing_gamma_prev;
std::vector<double> wake_gamma;
std::vector<Eigen::Vector3d> cpos, normal, tangent;

void InitParam() {
  AR = FLAGS_AR;
  ROWS = 4;
  COLS = AR * 3;
  CHORD = 1;
  SPAN = CHORD * AR;
  dx = CHORD / ROWS;
  dy = SPAN / COLS;
  INF = 1e20;
  alpha = FLAGS_alpha / 180. * M_PI;
  Q = FLAGS_Q;
  U = Eigen::Vector3d(Q, 0, 0);
  wing_gamma.resize((ROWS+1)*(COLS+1), 0);
}

template <class T>
auto& get_pos(T& v, std::size_t i, std::size_t j) {
  return v[j + i * (COLS + 1)];
}
template <class T>
const auto& get_pos(const T& v, std::size_t i, std::size_t j) {
  return v[j + i * (COLS + 1)];
}

std::size_t panel_index(std::size_t i,std::size_t j) {
  return j + i * COLS;
}
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
std::tuple<std::size_t, std::size_t> panel_index_inv(std::size_t index) {
  return std::make_tuple(index/COLS, index%COLS);
}

void InitPosition(std::vector<Eigen::Vector3d>& pos) {
  pos.resize((ROWS + 1) * (COLS + 1));
  // bound
  for (std::size_t i = 0; i < ROWS + 1; i++) {
    for (std::size_t j = 0; j < COLS + 1; j++) {
      double x0 = i * dx;
      double x = x0 * cos(-alpha);
      double y = -SPAN/2 + j * dy;
      double z = x0 * sin(-alpha);
      get_pos(pos, i, j) = Eigen::Vector3d(x, y, z);
    }
  }
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
      Eigen::Vector3d t =(
          (get_pos(pos, i + 1, j + 1) - get_pos(pos, i, j)) + 
          (get_pos(pos, i + 1, j) - get_pos(pos, i, j+1)))/2;
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
    const std::vector<double>& gammas,
    std::size_t i, std::size_t j) {
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
  for (std::size_t i=0; i<ROWS; ++i) {
    for (std::size_t j=0; j<COLS; ++j) {
      res += VORING(x, wing_pos, wing_gamma, i, j);
    }
  }
  return res;
}

Eigen::Vector3d WakeVelocity(const Eigen::Vector3d& x) {
  if (wake_gamma.size()==0) return Eigen::Vector3d::Zero();
  Eigen::Vector3d res = Eigen::Vector3d::Zero();
  std::size_t rows = wake_gamma.size() / COLS;
  for (std::size_t i=0; i<rows-1; ++i) {
    for (std::size_t j=0; j<COLS; ++j) {
      res += VORING(x, wake_pos, wake_gamma, i, j);
    }
  }
  return res;
}

auto CalcMatrix(const std::vector<Eigen::Vector3d>& cpos,
                const std::vector<Eigen::Vector3d>& normal) {
  // A_kl
  Eigen::MatrixXd res(ROWS * COLS, ROWS * COLS);
  std::fill(wing_gamma.begin(), wing_gamma.end(), 1);

  // loop for all vortices
  for (std::size_t K=0; K<ROWS*COLS; K++) {
    const auto cp = cpos[K];
    const auto n = normal[K];
    // loop for other vortices
    for (std::size_t L=0; L<ROWS*COLS; L++) {
      std::size_t i, j;
      std::tie(i, j) = panel_index_inv(L);
      auto u = VORING(cp, wing_pos, wing_gamma, i, j);
      res(K, L) = u.dot(n);
    }
  }
  return res;
}

auto CalcRhs(const std::vector<Eigen::Vector3d>& cpos,
             const std::vector<Eigen::Vector3d>& normal) {
  Eigen::VectorXd res(ROWS*COLS);
  std::size_t i, j;
  for (std::size_t K=0; K<ROWS*COLS; ++K) {
    std::tie(i, j) = panel_index_inv(K);
    Eigen::Vector3d u = U + WakeVelocity(cpos[K]);
    res(K) = -u.dot(normal[K]);
  }
  return res;
}

Eigen::Vector3d Velocity(const Eigen::Vector3d& x) {
  return U + BoundVelocity(x) + WakeVelocity(x);
}

Eigen::Vector3d CalcLift() {
  // TODO time derivation
  Eigen::Vector3d res = Eigen::Vector3d::Zero();
  for (std::size_t i=0; i<ROWS;i++) {
    for (std::size_t j=0; j<COLS;j++) {
      const auto cp = get_panel(cpos, i, j);
      Eigen::Vector3d Um = WakeVelocity(cp) + U;
      double g = get_panel(wing_gamma, i, j);
      double g_1 = i==0 ? 0 : get_panel(wing_gamma, i-1, j);
      double dy = fabs(get_pos(wing_pos, i, j).y() -get_pos(wing_pos, i, j+1).y());
      double dL = Um.dot(get_panel(tangent, i, j)) * (g -g_1) * dy;
      res += get_panel(normal, i, j) * dL;
    }
  }
  return res;
}

struct VortexLine {
  Eigen::Vector3d p0, p1;
  double g;
};

std::vector<VortexLine> GetLines() {
  std::vector<VortexLine> res;
  for (std::size_t i=0; i<ROWS;i++) {
    for (std::size_t j=0; j<COLS;j++) {
      std::vector<Eigen::Vector3d> corner = {
          get_pos(wing_pos, i, j), get_pos(wing_pos, i, j + 1),
          get_pos(wing_pos, i + 1, j + 1), get_pos(wing_pos, i + 1, j)};
      for (std::size_t k=0; k<corner.size(); k++) {
        if (i==ROWS-1 && k==2) continue;  // skip T.E
        res.push_back(VortexLine{corner[k], corner[(k + 1) % corner.size()],
                                 get_panel(wing_gamma, i, j)});
      }
    }
  }
  return res;
}

// Simpson's method
Eigen::Vector3d CalcLift2() {
  Eigen::Vector3d res = Eigen::Vector3d::Zero();
  // steady part
  auto lines = GetLines();
  for (const auto& line : lines) {
    Eigen::Vector3d mp = (line.p0 + line.p1) / 2;
    Eigen::Vector3d u = Velocity(mp);
    Eigen::Vector3d df = u.cross(line.p1 - line.p0) * line.g;
    res += df;
  }
  return res;
}

Eigen::Vector3d CalcLift2_unst(double dt) {
  // unsteady part
  Eigen::Vector3d res = Eigen::Vector3d::Zero();
  for (std::size_t i=0; i<ROWS;i++) {
    for (std::size_t j=0; j<COLS;j++) {
      const double b =
          (get_panel(wing_pos, i, j) - get_panel(wing_pos, i + 1, j)).norm();
      const double c =
          (get_panel(wing_pos, i, j) - get_panel(wing_pos, i, j + 1)).norm();
      const double A = b*c;
      const double dg_dt =
          (get_panel(wing_gamma, i, j) - get_panel(wing_gamma_prev, i, j)) / dt;
      LOG(INFO) << dg_dt << " " << A;
      res += dg_dt * A * get_panel(normal, i, j); 
    }
  }
  return res;
}

void MainLoop(std::size_t step, double dt) {
  wing_gamma_prev = wing_gamma;

  // TODO morphing
  
  // shed wake
  std::vector<Eigen::Vector3d> new_wake_pos;
  std::vector<double> new_wake_gamma;
  for (std::size_t j=0; j<=COLS; j++) {
    new_wake_pos.push_back(get_pos(wing_pos, ROWS, j) + U * dt);
  }
  if (step>1) {
    for (std::size_t j=0; j<COLS; j++) {
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
  auto A = CalcMatrix(cpos, normal);
  auto rhs = CalcRhs(cpos, normal);
  Eigen::FullPivLU<Eigen::MatrixXd> solver(A);
  Eigen::VectorXd gamma_v = solver.solve(rhs);
  for (std::size_t K = 0; K < ROWS * COLS; ++K) wing_gamma[K] = gamma_v(K);

  // calc load
  const auto L = CalcLift();
  LOG(INFO) << L.transpose();
  // const auto Fst = CalcLift2();
  // const auto Funst = CalcLift2_unst(dt);
  // LOG(INFO) << Fst.transpose() << " " << Funst.transpose();
  // const auto L = Fst + Funst;
  std::cout << L.z() / (0.5 * Q * Q * CHORD * SPAN) << std::endl;

  // output

  // advection
  std::vector<Eigen::Vector3d> wake_vel(wake_pos.size());
  std::transform(wake_pos.begin(), wake_pos.end(), wake_vel.begin(),
                 [](const auto& x) { return Velocity(x); });
  for (std::size_t i=0; i<wake_pos.size(); i++) {
    wake_pos[i] += wake_vel[i] * dt;
  }

  std::ofstream ofs("unsteady_problem.dat");
  CHECK(ofs) << "unable to open";
  // dump
  // 0
  for (auto p : wing_pos) ofs << p.transpose() << std::endl;
  ofs << std::endl << std::endl;
  // 1
  for (auto p : cpos) ofs << p.transpose() << std::endl;
  ofs << std::endl << std::endl;
  // 2
  ofs << "#normal";
  for (std::size_t i=0;i<ROWS*COLS;i++) {
    ofs << cpos[i].transpose() << "\t" << normal[i].transpose() << std::endl;
  }
  ofs << std::endl << std::endl;
  // 3 velocity at collocation points
  ofs << "#velocity at cp\n";
  for (std::size_t i=0;i<ROWS*COLS;i++) {
    ofs << cpos[i].transpose() << "\t" << Velocity(cpos[i]).transpose() << std::endl;
  }
  ofs << std::endl << std::endl;
  // 4 gamma
  ofs << "#gamma\n";
  for (std::size_t i=0; i<ROWS; ++i) {
    for (std::size_t j=0; j<COLS; ++j) {
      auto data = get_panel(cpos, i, j);
      data.z() = get_panel(wing_gamma, i, j);
      ofs << data.transpose() << std::endl;
    }
    ofs << std::endl;
  }
  ofs << std::endl << std::endl;
  // 5 rhs
  for (std::size_t i=0;i<ROWS*COLS;i++) {
    ofs << cpos[i].transpose() << "\t" << rhs[i] << std::endl;
  }
  ofs << std::endl << std::endl;
  // 6
  // self vorint gamma=1
  for (std::size_t i=0; i<ROWS; ++i) {
    for (std::size_t j=0; j<COLS; ++j) {
      auto cp = get_panel(cpos, i, j);
      auto u = Eigen::Vector3d::Zero();
      ofs << cp.transpose() << "\t" << u.transpose() << std::endl;
    }
  }
  ofs << std::endl << std::endl;
  // 7
  // wake velocity
  for (std::size_t i=0; i<=13; ++i) {
    for (std::size_t j=0; j<=13; ++j) {
      double x = 1.5;
      double ymin = -SPAN/2 - 1.;
      double ymax =  SPAN/2 + 1.;
      double zmin = -1.5;
      double zmax = 1.5;
      double y = ymin + (ymax-ymin) / 13 * i;
      double z = zmin + (zmax-zmin) / 13 * j;
      Eigen::Vector3d p(x, y, z);
      Eigen::Vector3d u = Velocity(p) - U;
      ofs << p.transpose() << "\t" << u.transpose() << std::endl;
    }
  }
  ofs << std::endl << std::endl;

  // 8 wake pos
  for (auto p : wake_pos) {
    ofs << p.transpose() << std::endl;
  }
  ofs << std::endl << std::endl;

  // 9 tangent
  ofs << "#normal";
  for (std::size_t i=0;i<ROWS*COLS;i++) {
    ofs << cpos[i].transpose() << "\t" << tangent[i].transpose() << std::endl;
  }
  ofs << std::endl << std::endl;
}

void SimulatorBody() { 
  InitPosition(wing_pos);
  double dt = 1./16.;
  for (std::size_t i=1; i<=50; i++) {
    LOG(INFO) << "step=" << i;
    MainLoop(i, dt);
  }
}

}  // anonymous namespace

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InstallFailureSignalHandler();
  FLAGS_logtostderr = true;
  InitParam();
  SimulatorBody();
  return 0;
}

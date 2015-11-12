/**
 * @file stady_problem_main.cpp
 * @brief Add description here
 *
 * output:
 *   stdout: C_L
 */

#include "../uvlm.h"
#include <fstream>
#include <gflags/gflags.h>
#include <glog/logging.h>

DEFINE_double(AR, 4, "aspect ratio");
DEFINE_double(Q, 1, "freestream velocity");
DEFINE_double(alpha, 5, "angle of attack [deg]");
DEFINE_bool(use_lift2, false, "calc lift in Simpson's way");

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

std::vector<Eigen::Vector3d> pos;
std::vector<double> gamma;

void InitParam() {
  AR = FLAGS_AR;
  ROWS = 4;
  COLS = AR * 4;
  CHORD = 1;
  SPAN = CHORD * AR;
  dx = CHORD / ROWS;
  dy = SPAN / COLS;
  INF = 1e20;
  alpha = FLAGS_alpha / 180. * M_PI;
  Q = FLAGS_Q;
  U = Eigen::Vector3d(Q, 0, 0);
  gamma.resize(ROWS*COLS);
}

template <class T>
auto& get_pos(T& v, std::size_t i, std::size_t j) {
  CHECK(v.size() == (COLS+1)*(ROWS+2)); return v[j + i * (COLS + 1)];
}
template <class T>
const auto& get_pos(const T& v, std::size_t i, std::size_t j) {
  CHECK(v.size() == (COLS+1)*(ROWS+2));
  return v[j + i * (COLS + 1)];
}

std::size_t panel_index(std::size_t i,std::size_t j) {
  return j + i * COLS;
}
template <class T>
auto& get_panel(T& v, std::size_t i, std::size_t j) {
  CHECK(v.size() == (COLS)*(ROWS));
  return v[j + i * COLS];
}
template <class T>
const auto& get_panel(const T& v, std::size_t i, std::size_t j) {
  CHECK(v.size() == (COLS)*(ROWS));
  return v[j + i * COLS];
}

void InitPosition() {
  pos.resize((ROWS + 2) * (COLS + 1));
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
  // wake
  for (std::size_t j = 0; j < COLS + 1; j++) {
    get_pos(pos, ROWS + 1, j) = Eigen::Vector3d(INF, j * dy, 0);
  }
}

auto CollocationPoints() {
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

auto Normals() {
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

auto VORTEX(const Eigen::Vector3d& x, const Eigen::Vector3d& x1,
    const Eigen::Vector3d& x2, double gamma) {
  // (10.16)
  // impl p. 584
  Eigen::Vector3d res;
  UVLM::BiotSavartLaw(&res, x1, x2, x);
  res *= gamma;
  return res;
}

Eigen::Vector3d VORING(const Eigen::Vector3d& x, std::size_t i,std::size_t j,
    double gamma) {
  Eigen::Vector3d u = Eigen::Vector3d::Zero();
  u += VORTEX(x, get_pos(pos, i, j), get_pos(pos, i, j + 1), gamma);
  u += VORTEX(x, get_pos(pos, i, j + 1), get_pos(pos, i + 1, j + 1), gamma);
  u += VORTEX(x, get_pos(pos, i + 1, j + 1), get_pos(pos, i + 1, j), gamma);
  u += VORTEX(x, get_pos(pos, i + 1, j), get_pos(pos, i, j), gamma);
  return u;
}

auto CalcMatrix(const std::vector<Eigen::Vector3d>& cpos,
                const std::vector<Eigen::Vector3d>& normal) {
  // A_kl
  Eigen::MatrixXd res(ROWS * COLS, ROWS * COLS);

  // loop for all bound vortices
  for (std::size_t i=0; i<ROWS; ++i) {
    for (std::size_t j=0; j<COLS; ++j) {
      const auto k = panel_index(i, j);
      const auto cp = get_panel(cpos, i, j);
      const auto n = get_panel(normal, i, j);
      // std::cerr << i << " " << j << " " << k << std::endl;
      
      // loop for vortex ring
      for (std::size_t ii=0; ii<ROWS; ++ii) {
        for (std::size_t jj=0; jj<COLS; ++jj) {
          const auto l = panel_index(ii, jj);
          auto u = VORING(cp, ii, jj, 1);
          res(k, l) = u.dot(n);
        }
      }
      // add influence of wake
      for (std::size_t jj=0; jj<COLS; ++jj) {
        const std::size_t ii = ROWS-1;
        const auto l = panel_index(ii, jj);
        auto u = VORING(cp, ii+1, jj, 1);
        res(k, l) += u.dot(n);
      }
    }
  }
  return res;
}

auto CalcRhs(const std::vector<Eigen::Vector3d>& normal) {
  Eigen::VectorXd res(ROWS*COLS);
  for (std::size_t i=0; i<ROWS*COLS; ++i) {
    res(i) = - U.dot(normal[i]);
  }
  return res;
}

Eigen::Vector3d Velocity(const Eigen::Vector3d& x) {
  Eigen::Vector3d res = U;
  // bound
  for (std::size_t i=0; i<ROWS; ++i) {
    for (std::size_t j=0; j<COLS; ++j) {
      res += VORING(x, i, j, get_panel(gamma, i, j));
    }
  }
  // wake
  for (std::size_t j=0; j<COLS; ++j) {
    std::size_t i = ROWS - 1;
    res += VORING(x, i+1, j, get_panel(gamma, i, j));
  }
  return res;
}

double CalcLift() {
  double res = 0;
  for (std::size_t i=0; i<ROWS;i++) {
    for (std::size_t j=0; j<COLS;j++) {
      double g = get_panel(gamma, i, j);
      double g_1 = i==0 ? 0 : get_panel(gamma, i-1, j);
      double dy = fabs(get_pos(pos, i, j).y() -get_pos(pos, i, j+1).y());
      double dL = Q * (g -g_1) * dy;
      res += dL;
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
          get_pos(pos, i, j), get_pos(pos, i, j + 1),
          get_pos(pos, i + 1, j + 1), get_pos(pos, i + 1, j)};
      for (std::size_t k=0; k<corner.size(); k++) {
        if (i==ROWS-1 && k==2) continue;  // skip T.E
        res.push_back(VortexLine{corner[k], corner[(k + 1) % corner.size()],
                                 get_panel(gamma, i, j)});
      }
    }
  }
  // add line in wake
  // std::size_t i = ROWS-1;
  // for (std::size_t j = 0; j < COLS; j++) {
  //   res.push_back(VortexLine{get_pos(pos, i + 1, j), get_pos(pos, i + 1, j + 1),
  //                            get_panel(gamma, i, j)});
  // }
  return res;
}

// Simpson's method
Eigen::Vector3d CalcLift2() {
  Eigen::Vector3d res = Eigen::Vector3d::Zero();
  auto lines = GetLines();
  for (auto line : lines) {
    Eigen::Vector3d mp = (line.p0 + line.p1) / 2;
    Eigen::Vector3d u = Velocity(mp);
    Eigen::Vector3d df = u.cross(line.p1 - line.p0) * line.g;
    res += df;
  }
  return res;
}

void SimulatorBody() { 
  InitPosition();
  const auto cpos = CollocationPoints();
  const auto normal = Normals();

  auto A = CalcMatrix(cpos, normal);
  auto rhs = CalcRhs(normal);
  Eigen::FullPivLU<Eigen::MatrixXd> solver(A);
  Eigen::VectorXd gamma_v = solver.solve(rhs);

  for (std::size_t i=0; i<ROWS*COLS;++i) gamma[i] = gamma_v(i);

  std::vector<Eigen::Vector3d> v;
  for (const auto& cp : cpos) 
    v.emplace_back(Velocity(cp));

  std::ofstream ofs("steady_problem.dat");
  CHECK(ofs) << "unable to open";
  // dump
  // 0
  for (auto p : pos) ofs << p.transpose() << std::endl;
  ofs << std::endl << std::endl;
  // 1
  for (auto p : cpos) ofs << p.transpose() << std::endl;
  ofs << std::endl << std::endl;
  // 2
  ofs << "normal";
  for (std::size_t i=0;i<ROWS*COLS;i++) {
    ofs << cpos[i].transpose() << "\t" << normal[i].transpose() << std::endl;
  }
  ofs << std::endl << std::endl;
  // 3 velocity at collocation points
  ofs << "velocity\n";
  for (std::size_t i=0;i<ROWS*COLS;i++) {
    ofs << cpos[i].transpose() << "\t" << v[i].transpose() << std::endl;
  }
  ofs << std::endl << std::endl;
  // 4 gamma
  ofs << "gamma\n";
  for (std::size_t i=0; i<ROWS; ++i) {
    for (std::size_t j=0; j<COLS; ++j) {
      auto data = get_panel(cpos, i, j);
      data.z() = get_panel(gamma, i, j);
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
      auto u = VORING(cp, i, j, 1);
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

  // 8 vortex lines
  auto lines = GetLines();
  for (auto line : lines) {
    ofs << line.p0.transpose() << std::endl;
  }

  // 9 mid point
  for (auto line : lines) {
    Eigen::Vector3d mp = (line.p0 + line.p1) / 2;
    Eigen::Vector3d u = Velocity(mp);
    ofs << mp.transpose() << "\t" << u.transpose() << std::endl;
  }

  // std::cout << "Lift=" << CalcLift() << std::endl;
  std::cout << CalcLift() / (0.5 * Q * Q * CHORD * SPAN) << std::endl;
  auto C = CalcLift2() / (0.5 * Q * Q * CHORD * SPAN);
  LOG(INFO) << C.x() << " " << C.z();
  //
  // double CL_ans = 2. * M_PI * sin(alpha);
  // std::cout << "CL(2D)=" << CL_ans << std::endl;
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

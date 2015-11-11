/**
 * @file stady_problem_main.cpp
 * @brief Add description here
 */

#include "../uvlm.h"

namespace {

const std::size_t ROWS = 4;
const std::size_t COLS = 26;
const double CHORD = 1;
const double AR = 6; 
const double SPAN = CHORD * AR;
const double dx = CHORD / ROWS;
const double dy = SPAN / COLS;
const double INF = 1e20;
const double alpha = 5. / 180. * M_PI;
const double Q = 1;
const Eigen::Vector3d U(Q, 0, 0);

std::vector<Eigen::Vector3d> pos;
std::vector<double> gamma(ROWS* COLS);

template <class T>
auto& get_pos(T& v, std::size_t i, std::size_t j) {
  CHECK(v.size() == (COLS+1)*(ROWS+2));
  return v[j + i * (COLS + 1)];
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
      double x = i * dx * cos(-alpha);
      double y = j * dy;
      double z = x * sin(-alpha);
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

  // dump
  // 0
  for (auto p : pos) std::cout << p.transpose() << std::endl;
  std::cout << std::endl << std::endl;
  // 1
  for (auto p : cpos) std::cout << p.transpose() << std::endl;
  std::cout << std::endl << std::endl;
  // 2
  for (std::size_t i=0;i<ROWS*COLS;i++) {
    std::cout << cpos[i].transpose() << "\t" << normal[i].transpose() << std::endl;
  }
  std::cout << std::endl << std::endl;
  // 3 velocity at collocation points
  std::cout << "velocity\n";
  for (std::size_t i=0;i<ROWS*COLS;i++) {
    std::cout << cpos[i].transpose() << "\t" << v[i].transpose() << std::endl;
  }
  std::cout << std::endl << std::endl;
  // 4 gamma
  std::cout << "gamma\n";
  for (std::size_t i=0; i<ROWS; ++i) {
    for (std::size_t j=0; j<COLS; ++j) {
      auto data = get_panel(cpos, i, j);
      data.z() = get_panel(gamma, i, j);
      std::cout << data.transpose() << std::endl;
    }
    std::cout << std::endl;
  }
  std::cout << std::endl << std::endl;
  // 5 rhs
  for (std::size_t i=0;i<ROWS*COLS;i++) {
    std::cout << cpos[i].transpose() << "\t" << rhs[i] << std::endl;
  }
  std::cout << std::endl << std::endl;
  // 6
  // self vorint gamma=1
  for (std::size_t i=0; i<ROWS; ++i) {
    for (std::size_t j=0; j<COLS; ++j) {
      auto cp = get_panel(cpos, i, j);
      auto u = VORING(cp, i, j, 1);
      std::cout << cp.transpose() << "\t" << u.transpose() << std::endl;
    }
  }
  std::cout << std::endl << std::endl;

  std::cerr << "Lift=" << CalcLift() << std::endl;
  std::cerr << "CL=" << CalcLift() / (0.5 * Q * Q * CHORD * SPAN) << std::endl;

  double CL_ans = 2. * M_PI * sin(alpha);
  std::cerr << "CL(2D)=" << CL_ans << std::endl;
}


}  // anonymous namespace

int main(int argc, char* argv[]) {
  SimulatorBody();
  return 0;
}

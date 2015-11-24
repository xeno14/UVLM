/**
 * @file stady_problem2_main.cpp
 * @brief Add description here
 */

#include "../uvlm.h"

#include <gflags/gflags.h>
#include <glog/logging.h>
#include <fstream>

DEFINE_double(AR, 4, "aspect ratio");
DEFINE_double(Q, 1, "freestream velocity");
DEFINE_double(alpha, 5, "angle of attack [deg]");

using namespace UVLM;

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

auto vortices = std::make_shared<std::vector<UVLM::VortexRing>>();
std::vector<UVLM::VortexContainer> containers;
std::vector<UVLM::Morphing> morphings;

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
}

auto InitWing() {
  UVLM::proto::Wing wing;
  UVLM::wing::RectGenerator()
      .Generate(&wing, CHORD, SPAN / 2, ROWS, COLS/2);
  UVLM::wing::SetOrigin(&wing, {0, 0, 0});
  UVLM::WingBuilder builder(&containers, vortices);
  builder.AddWing(wing);
  builder.Build();
}

void AddWake() {
  std::vector<VortexRing> wake(containers[0].edge_begin(),
                               containers[0].edge_end());
  for (auto& w : wake) {
    auto& node = w.nodes();
    node[0] = node[1];
    node[3] = node[2];
    node[1].x() = INF;
    node[2].x() = INF;
  }
  vortices->insert(vortices->end(), wake.begin(), wake.end());
}

void SimulatorBody() {
  const std::size_t wake_offset = ROWS*COLS;
  InitWing();
  Morphing m;
  m.set_alpha(alpha);
  auto& c = containers[0];
  for (auto& v : c) {
    for (std::size_t i = 0; i < v.nodes().size(); i++) {
      m.Perfome(&v.nodes()[i], v.nodes0()[i], 0);
    }
  }
  AddWake();

  morphings.push_back(m);

  int n=10;
  while(n--) {
    auto gamma = SolveLinearProblem(containers, morphings, U, 0);
    for (std::size_t i = 0; i < wake_offset; i++) {
      vortices->at(i).set_gamma(gamma(i));
    }
    std::size_t i=0;
    for (auto it=c.edge_begin(); it!=c.edge_end(); ++it) {
      (vortices->begin() + c.size() + i)->set_gamma(it->gamma());
    }
  }
}

Eigen::Vector3d Velocity(const Eigen::Vector3d& p) {
  Eigen::Vector3d u;
  InducedVelocity(&u, p, vortices->begin(), vortices->end());
  return u + U;
}

}  // anonymous namespace

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InstallFailureSignalHandler();
  FLAGS_logtostderr = true;
  InitParam();
  SimulatorBody();

  std::ofstream ofs("steady_problem2.dat");
  CHECK(ofs) << "unable to open";
  // 0 points
  for (auto v : *vortices) {
    for (auto n : v.nodes()) ofs<<n.transpose() <<std::endl;
  }
  ofs << std::endl << std::endl;

  const auto& c = containers[0];

  // 1 velocity at centroid
  for (auto v : c) {
    auto cp = v.Centroid();
    auto u = Velocity(cp);
    ofs << cp.transpose() << "\t" << u.transpose() << std::endl;
  }
  ofs << std::endl << std::endl;

  // 2
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

  // 3 gamma
  for (auto v : c) {
    auto cp = v.Centroid();
    cp.z() = v.gamma();
    ofs << cp.transpose()<<std::endl;
  }
  ofs << std::endl << std::endl;

  return 0;
}

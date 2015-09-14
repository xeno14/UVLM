/**
 * @file plug_main.cpp
 * @brief Add description here
 */


#include "../proto/uvlm.pb.h"
// #include "calc_load.h"
#include "uvlm_vortex_ring.h"
#include "morphing.h"
#include "linear.h"
#include "shed.h"
#include "proto_adaptor.h"
#include "util.h"
#include "vortex_container.h"
#include "wing/wing.h"
#include "wing_builder.h"
#include "parameter.h"

#include <gflags/gflags.h>
#include <glog/logging.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <memory>

DEFINE_string(input, "", "setting yaml");
DEFINE_string(output, "", "output path");

namespace {

YAML::Node config;

}  // anonymous namespace

std::vector<Eigen::Vector3d> ReadWingTsv (const std::string& path) {
  std::vector<Eigen::Vector3d> res;
  std::ifstream ifs(path);
  if (!ifs) {
    LOG(INFO) << "Open error" ;
    std::exit(EXIT_FAILURE);
  }
  while (!ifs.eof()) {
    double x, y, z;
    ifs >> x >> y >> z;
    res.emplace_back(x, y, z);
  }
  return res;
}

auto InitWing(UVLM::WingBuilder* builder,
              const std::vector<Eigen::Vector3d>& origins) {
  UVLM::proto::Wing wing;

  // 翼の作成
  // NACA0012: AR=8
  DEFINE_PARAM(int, rows, config["parameter"]);
  DEFINE_PARAM(int, cols, config["parameter"]);
  UVLM::wing::NACA00XXGenerator(12, 1., 4., PARAM_rows, PARAM_cols)
      .Generate(&wing);

  for (const auto& origin : origins) {
    auto* p = wing.mutable_origin();
    p->CopyFrom(UVLM::Vector3dToPoint(origin));
    builder->AddWing(wing);
  }
  builder->Build();

  return wing;
}

template <class InputIterator1, class InputIterator2>
void OutputSnapshot2(const std::size_t index, InputIterator1 container_first,
                    InputIterator1 container_last, InputIterator2 vortex_first,
                    InputIterator2 vortex_last, const double t) {
  UVLM::proto::Snapshot2 snapshot;
  snapshot.set_t(t);

  for(auto container = container_first; container != container_last; ++container) {
    auto* shape = snapshot.add_container_shapes();
    shape->set_rows(container->rows());
    shape->set_cols(container->cols());
    shape->set_id(container->id());
    // shape->set_origin(Vector3dToPoint(container->origin()));
  }
  for (auto vortex = vortex_first; vortex != vortex_last; ++vortex) {
    snapshot.add_vortices()->CopyFrom(VortexRingToProto(*vortex));
  }
  char filename[256];
  sprintf(filename, "%s/%08lu", FLAGS_output.c_str(), index); 
  std::ofstream ofs(filename, std::ios::binary);
  CHECK_OPEN(ofs);
  snapshot.SerializeToOstream(&ofs);
}


void SimulationBody() {
  auto vortices = std::make_shared<std::vector<UVLM::VortexRing>>();
  UVLM::UVLMVortexRing rings;
  UVLM::Morphing morphing;

  const auto param = config["parameter"];
  DEFINE_PARAM(double, alpha, param);
  DEFINE_PARAM(double, U, param);
  DEFINE_PARAM(double, dt, param);

  const double ALPHA = M_PI * 2 / 360 * PARAM_alpha;  // Angle of attack
  const double U = PARAM_U;                           // Upstream velocity
  const double K = 0.1;                               // Reduced frequency
  const double C = 1;                                 // Chord length
  const double OMEGA = 2 * U * K / C;                 // Flapping frequency
  const double PHI = 15 * M_PI / 180;                 // Angle of flapping
  // const double T = M_PI * 2 / OMEGA;                  // Period of flapping
  const double BETA = 4 * M_PI / 180;                 // twising amp at wing tip
  // const double RHO = 1;                               // Fluid density
  Eigen::Vector3d Vinfty(U * cos(ALPHA), 0, U * sin(ALPHA));

  std::vector<UVLM::VortexContainer> containers;
  std::vector<UVLM::VortexContainer> containers_prev;
  std::vector<Eigen::Vector3d> origins { {0, 0, 0} };
  std::vector<UVLM::Morphing> morphings;

  UVLM::WingBuilder builder(&containers, vortices);
  auto wing = InitWing(&builder, origins);

  const double CHORD = wing.chord();    // chord length
  const double SPAN = wing.span();      // wing span
  LOG(INFO) << "chord: " << CHORD;
  LOG(INFO) << "span: " << SPAN;

  containers_prev.resize(containers.size());

  rings.bound_vortices() = *vortices;

  // morphing.set_plug([OMEGA](double t) { return 0.2 * sin(OMEGA*t); });
  morphing.set_flap([OMEGA, PHI](double t) { return PHI * sin(OMEGA * t + M_PI); });
  morphing.set_twist([OMEGA, BETA, SPAN](const Eigen::Vector3d& x0, double t) {
    return BETA * fabs(x0.y()) / SPAN * sin(OMEGA * t + M_PI);
  });

  const double dt = PARAM_dt;

  std::size_t wake_offset =
      CountTotalSize(containers.cbegin(), containers.cend());

  DEFINE_PARAM(int, steps, config["setting"]);
  // main loop
  for(int ti=0; ti<PARAM_steps; ti++) {
    LOG(INFO) << "step: " << ti;
    const double t = ti * dt;

    // 連立方程式を解いて翼の上の循環を求める
    LOG(INFO) << "Solve linear problem" ;
    auto gamma = UVLM::SolveLinearProblem(
        vortices->begin(), vortices->begin() + wake_offset,
        vortices->cbegin() + wake_offset, vortices->cend(),
        Vinfty, morphing, t);
    for (std::size_t i=0; i<wake_offset; i++) {
      vortices->at(i).set_gamma(gamma(i));
    }

    // TODO remove rings
    LOG(INFO) << "Output" ;
    std::copy(vortices->begin(), vortices->begin() + wake_offset,
              rings.bound_vortices().begin());
    OutputSnapshot2(ti, containers.begin(), containers.end(), vortices->begin(),
                   vortices->end(), t);

    LOG(INFO) << "Shed" ;
    // 放出する渦を求める
    // TODO 変形したときに位置が合わない
    // TODO remove rings
    std::vector<std::vector<UVLM::VortexRing>> shed(containers.size());

    // TODO use zip iterator?
    for (std::size_t i=0; i<containers.size(); i++) {
      auto edge_first = containers[i].edge_begin();
      auto edge_last = containers[i].edge_end();
      shed[i].resize(std::distance(edge_first, edge_last));
      UVLM::ShedAtTrailingEdge(edge_first, edge_last, shed[i].begin(),
                               vortices->cbegin(), vortices->cend(), 
                               rings,
                               Vinfty, t,
                               dt);
    }

    LOG(INFO) << "Advect" ;
    // TODO remove rings
    UVLM::AdvectWake(vortices->begin() + wake_offset, vortices->end(),
                     vortices->cbegin(), vortices->cend(), 
                     rings,
                     Vinfty, dt);

    LOG(INFO) << "Morphing" ;
    // 変形する
    for (std::size_t i=0; i<containers.size(); i++) {
      for (auto& vortex : containers[i]) {
        for (std::size_t j=0; j<vortex.nodes().size(); j++) {
          morphing.Perfome(&vortex.nodes()[j], vortex.nodes0()[j], t, dt);
        }
      }
    }

    LOG(INFO) << "Edge" ;
    // 変形後のedgeの位置とshedを合わせる
    for (std::size_t i=0; i<containers.size(); i++) {
      UVLM::AttachShedVorticesToEdge(containers[i].edge_begin(), containers[i].edge_end(),
                                     shed[i].begin());
    }

    LOG(INFO) << "Append shed" ;
    // 放出した渦の追加
    // TODO jointed iterator?
    for (std::size_t i=0; i<shed.size(); i++) {
      vortices->insert(vortices->end(), shed[i].cbegin(), shed[i].cend());
    }

    LOG(INFO) << "copy" ;
    rings.wake_vortices().resize(vortices->size() - wake_offset);
    std::copy(vortices->begin() + wake_offset, vortices->end(),
              rings.wake_vortices().begin());

    std::cerr << std::endl;

    // Calc aerodynamic loads
    // LOG(INFO) << "aerodynamic load";
    // if (ti >= 1) {
    //   for (std::size_t i = 0; i < containers.size(); ++i) {
    //     auto force = UVLM::CalcLoad(containers[i], containers_prev[i],
    //                                 vortices->begin() + wake_offset,
    //                                 vortices->end(), Vinfty, RHO, dt);
    //     Eigen::Vector3d coeff = force / (0.5 * RHO * U * U);
    //     printf("%e\t%e\t%e\t%e\n", t/T, coeff.x(), coeff.y(), coeff.z());
    //   }
    // }
    // LOG(INFO) << "Copy containers";
    // CopyContainers(containers.cbegin(), containers.cend(),
    //                containers_prev.begin());
  }
}

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  FLAGS_logtostderr = true;   // TODO 実行時に --logtostderr にすると怒られる

  config = YAML::LoadFile(FLAGS_input);
#ifdef _OPENMP
  auto setting = config["setting"];
  DEFINE_PARAM(int, threads, setting);
  LOG(INFO) << "work with " << PARAM_threads << " threads." ;
  omp_set_num_threads(PARAM_threads);
#endif
  SimulationBody();
  return 0;
}

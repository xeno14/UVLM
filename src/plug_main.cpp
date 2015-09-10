/**
 * @file plug_main.cpp
 * @brief Add description here
 */


#include "../proto/uvlm.pb.h"
#include "calc_load.h"
#include "uvlm_vortex_ring.h"
#include "morphing.h"
#include "linear.h"
#include "shed.h"
#include "proto_adaptor.h"
#include "util.h"
#include "vortex_container.h"
#include "wing_builder.h"

#include <gflags/gflags.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <memory>

DEFINE_string(output, "", "output path");
DEFINE_string(wing, "", "wing data");
DEFINE_double(dt, 0.01, "delta t");
#ifdef _OPENMP
DEFINE_int32(threads, 1, "openmp threads");
#endif
DEFINE_int32(steps, 100, "number of time steps");

std::vector<Eigen::Vector3d> ReadWingTsv (const std::string& path) {
  std::vector<Eigen::Vector3d> res;
  std::ifstream ifs(path);
  if (!ifs) {
    std::cerr << "Open error" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  while (!ifs.eof()) {
    double x, y, z;
    ifs >> x >> y >> z;
    res.emplace_back(x, y, z);
  }
  return res;
}

void InitWing(UVLM::WingBuilder* builder,
              const std::vector<Eigen::Vector3d>& origins) {
  std::ifstream ifs(FLAGS_wing, std::ios::binary);
  CHECK_OPEN(ifs);
  UVLM::proto::Wing wing;
  wing.ParseFromIstream(&ifs);
  for (const auto& origin : origins) {
    auto* p = wing.mutable_origin();
    p->CopyFrom(UVLM::Vector3dToPoint(origin));
    builder->AddWing(wing);
  }
  builder->Build();
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

std::size_t CountTotalSize(const std::vector<UVLM::VortexContainer>& containers) {
  std::size_t res = 0;
  for (const auto& c : containers) res += c.size();
  return res;
}

void CopyPrevContainers(std::vector<UVLM::VortexContainer>* const prev,
                      const std::vector<UVLM::VortexContainer>& containers) {
  auto vortices = containers.begin()->vortices();
  const auto total_size = CountTotalSize(containers);
  auto vortices_prev = std::make_shared<std::vector<UVLM::VortexRing>>(
      vortices->begin(), vortices->begin() + total_size);
  prev->clear();
  for (const auto& c : containers) {
    prev->emplace_back(vortices_prev, c.rows(), c.cols(), c.id());
  }
}

void SimulationBody() {
  auto vortices = std::make_shared<std::vector<UVLM::VortexRing>>();
  UVLM::UVLMVortexRing rings;
  UVLM::Morphing morphing;
  const double ALPHA = M_PI * 2 / 360 * 4;  // Angle of attack
  const double U = 1;                       // Upstream velocity
  const double K = 0.1;                     // Reduced frequency
  const double C = 1;                       // Chord length
  const double OMEGA = 2 * U * K / C;       // Flapping frequency
  const double PHI = M_PI * 2 / 360 * 15;   // Angle of flapping
  Eigen::Vector3d Vinfty(U * cos(ALPHA), 0, sin(ALPHA));

  std::vector<UVLM::VortexContainer> containers;
  std::vector<UVLM::VortexContainer> containers_prev;
  std::vector<Eigen::Vector3d> origins;
  std::vector<UVLM::Morphing> morphings;

  origins.emplace_back(0, 0, 0);

  UVLM::WingBuilder builder(&containers, vortices);
  InitWing(&builder, origins);

  rings.bound_vortices() = *vortices;

  // morphing.set_plug([OMEGA](double t) { return 0.2 * sin(OMEGA*t); });
  morphing.set_flap([OMEGA, PHI](double t) { return PHI * sin(OMEGA * t); });

  const double dt = FLAGS_dt;

  std::size_t wake_offset = CountTotalSize(containers);

  std::cerr << vortices->size() <<"aa\n";
  // main loop
  for(int ti=0; ti<FLAGS_steps; ti++) {
    std::cerr << ti << std::endl;
    const double t = ti * dt;

    // 連立方程式を解いて翼の上の循環を求める
    std::cerr << "Solve linear problem" << std::endl;
    auto gamma = UVLM::SolveLinearProblem(
        vortices->begin(), vortices->begin() + wake_offset,
        vortices->cbegin() + wake_offset, vortices->cend(),
        Vinfty, morphing, t);
    for (std::size_t i=0; i<wake_offset; i++) {
      vortices->at(i).set_gamma(gamma(i));
    }

    // TODO remove rings
    std::cerr << "Output" << std::endl;
    std::copy(vortices->begin(), vortices->begin() + wake_offset,
              rings.bound_vortices().begin());
    OutputSnapshot2(ti, containers.begin(), containers.end(), vortices->begin(),
                   vortices->end(), t);

    std::cerr << "Shed" << std::endl;
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

    std::cerr << "Advect" << std::endl;
    // TODO remove rings
    UVLM::AdvectWake(vortices->begin() + wake_offset, vortices->end(),
                     vortices->cbegin(), vortices->cend(), 
                     rings,
                     Vinfty, dt);

    std::cerr << "Morphing" << std::endl;
    // 変形する
    for (std::size_t i=0; i<containers.size(); i++) {
      for (auto& vortex : containers[i]) {
        for (std::size_t j=0; j<vortex.nodes().size(); j++) {
          morphing.Perfome(&vortex.nodes()[j], vortex.nodes0()[j], t, dt);
        }
      }
    }

    std::cerr << "Edge" << std::endl;
    // 変形後のedgeの位置とshedを合わせる
    for (std::size_t i=0; i<containers.size(); i++) {
      UVLM::AttachShedVorticesToEdge(containers[i].edge_begin(), containers[i].edge_end(),
                                     shed[i].begin());
    }

    std::cerr << "Append shed" << std::endl;
    // 放出した渦の追加
    // TODO jointed iterator?
    for (std::size_t i=0; i<shed.size(); i++) {
      vortices->insert(vortices->end(), shed[i].cbegin(), shed[i].cend());
    }

    std::cerr << "copy" << std::endl;
    // TODO remove rings
    // std::copy(advected_wake.begin(), advected_wake.end(),
    //           vortices->begin() + wake_offset);
    rings.wake_vortices().resize(vortices->size() - wake_offset);
    std::copy(vortices->begin() + wake_offset, vortices->end(),
              rings.wake_vortices().begin());
    // for (std::size_t i = 0; i < container.size(); i++) {
    //   std::cerr << vortices->at(i).gamma() << " vs " << container[i].gamma()
    //             << " vs " << rings.bound_vortices()[i].gamma() << std::endl;
    // }

    std::cerr << std::endl;

    CopyPrevContainers(&containers_prev, containers);
    if (ti >= 1) {
      for (std::size_t i = 0; i < containers.size(); ++i) {
        auto force = UVLM::CalcLoad(containers[i], containers_prev[i],
                                    vortices->begin() + wake_offset,
                                    vortices->end(), Vinfty, 1, dt);
        printf("%e\t%e\t%e\t%e\n", t, force.x(), force.y(), force.z());
      }
    }
  }
}

int main(int argc, char* argv[]) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
#ifdef _OPENMP
  std::cerr << "work with " << FLAGS_threads << " threads." << std::endl;
  omp_set_num_threads(FLAGS_threads);
#endif
  SimulationBody();
  return 0;
}

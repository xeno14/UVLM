/**
 * @file plug_main.cpp
 * @brief Add description here
 */


#include "../proto/uvlm.pb.h"
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

void InitWing(UVLM::UVLMVortexRing* rings) {
  std::ifstream ifs(FLAGS_wing, std::ios::binary);
  CHECK_OPEN(ifs);
  UVLM::proto::Wing wing;
  wing.ParseFromIstream(&ifs);
  std::vector<Eigen::Vector3d> pos(wing.points().size());
  std::transform(wing.points().begin(), wing.points().end(), pos.begin(),
                 UVLM::PointToVector3d);
  rings->InitWing(pos, wing.cols());
}

void InitWing(UVLM::WingBuilder* builder) {
  std::ifstream ifs(FLAGS_wing, std::ios::binary);
  CHECK_OPEN(ifs);
  UVLM::proto::Wing wing;
  wing.ParseFromIstream(&ifs);
  builder->AddWing(wing).Build();
}


void OutputSnapshot(const int index, const double t, const UVLM::UVLMVortexRing& rings) {
  UVLM::proto::Snapshot snapshot;
  snapshot.set_t(t);

  auto* flying_wing = snapshot.add_flying_wings();
  UVLMVortexRingToBird(flying_wing, rings);

  char filename[256];
  sprintf(filename, "%s/%08d", FLAGS_output.c_str(), index); 
  std::ofstream ofs(filename, std::ios::binary);
  if (!ofs) {
    std::cerr << "Open error " << filename << std::endl;
    std::exit(EXIT_FAILURE);
  }
  snapshot.SerializeToOstream(&ofs);
}

void SimulationBody() {
  auto vortices = std::make_shared<std::vector<UVLM::VortexRing>>();
  UVLM::UVLMVortexRing rings;
  UVLM::Morphing morphing;
  Eigen::Vector3d Vinfty(2, 0, 0.1);

  std::vector<UVLM::VortexContainer> containers;

  UVLM::WingBuilder builder(&containers, vortices);
  InitWing(&builder);

  rings.bound_vortices() = *vortices;

  auto& container = containers[0];

  morphing.set_plug([](double t) { return 0.2 * sin(5*t); });
  morphing.set_flap([](double t) { return M_PI/6 * sin(4*t); });

  const double dt = FLAGS_dt;

  const std::size_t wake_offset = container.size();

  std::cerr << vortices->size() <<"aa\n";
  // main loop
  for(std::size_t i=0; i<100; i++) {
    std::cerr << i << std::endl;
    const double t = i * dt;

    // 連立方程式を解いて翼の上の循環を求める
    std::cerr << "Solve linear problem" << std::endl;
    auto gamma = UVLM::SolveLinearProblem(
        container.begin(), container.end(),
        vortices->cbegin() + wake_offset, vortices->cend(),
        Vinfty, morphing, t);
    for (std::size_t i=0; i<wake_offset; i++) {
      vortices->at(i).set_gamma(gamma(i));
    }

    // TODO remove rings
    std::cerr << "Output" << std::endl;
    std::copy(container.begin(), container.end(),
              rings.bound_vortices().begin());
    OutputSnapshot(i, t, rings);

    std::cerr << "Shed" << std::endl;
    // 放出する渦を求める
    // TODO 変形したときに位置が合わない
    // TODO remove rings
    auto edge_first = container.edge_begin();
    auto edge_last = container.edge_end();
    std::vector<UVLM::VortexRing> shed(std::distance(edge_first, edge_last));
    UVLM::ShedAtTrailingEdge(edge_first, edge_last, shed.begin(),
                             vortices->cbegin(), vortices->cend(), 
                             rings,
                             Vinfty, t,
                             dt);

    std::cerr << "Advect" << std::endl;
    // TODO remove rings
    UVLM::AdvectWake(vortices->begin() + wake_offset, vortices->end(),
                     vortices->cbegin(), vortices->cend(), 
                     rings,
                     Vinfty, dt);

    std::cerr << "Morphing" << std::endl;
    // 変形する
    for (auto& vortex : container) {
      for (std::size_t i=0; i<vortex.nodes().size(); i++) {
        morphing.Perfome(&vortex.nodes()[i], vortex.nodes0()[i], t, dt);
      }
    }

    std::cerr << "Edge" << std::endl;
    // 変形後のedgeの位置とshedを合わせる
    UVLM::AttachShedVorticesToEdge(container.edge_begin(), container.edge_end(),
                                   shed.begin());

    std::cerr << "Append shed" << std::endl;
    // 放出した渦の追加
    vortices->insert(vortices->end(), shed.cbegin(), shed.cend());

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
  }
}

int main(int argc, char* argv[]) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  SimulationBody();
  return 0;
}

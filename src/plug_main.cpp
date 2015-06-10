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
#include "vortex_container.h"

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
  if (!ifs) {
    std::cerr << "Cannot open " << FLAGS_wing << std::endl;
    std::exit(EXIT_FAILURE);
  }
  UVLM::proto::Wing wing;
  wing.ParseFromIstream(&ifs);
  std::vector<Eigen::Vector3d> pos(wing.points().size());
  std::transform(wing.points().begin(), wing.points().end(), pos.begin(),
                 UVLM::PointToVector3d);
  rings->InitWing(pos, wing.cols());
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

  const Eigen::Vector3d origin(0, 0, 0);
  morphing.set_origin(origin);
  rings.set_origin(origin);

  // TODO ringsを取り除く
  InitWing(&rings);
  *vortices = rings.bound_vortices();

  UVLM::VortexContainer container(vortices, rings.rows(), rings.cols() * 2, 0);

  morphing.set_plug([](double t) { return 0.2 * sin(5*t); });
  // morphing.set_flap([](double t) { return M_PI/6 * sin(4*t); });

  const double dt = FLAGS_dt;

  const std::size_t wake_offset = container.size();

  // main loop
  for(std::size_t i=0; i<100; i++) {
    std::cerr << i << std::endl;
    const double t = i * dt;

    // 連立方程式を解いて翼の上の循環を求める
    auto gamma = UVLM::SolveLinearProblem(
        container.begin(), container.end(),
        vortices->cbegin() + wake_offset, vortices->cend(),
        Vinfty, morphing, t);
    for (std::size_t i=0; i<wake_offset; i++) {
      vortices->at(i).set_gamma(gamma(i));
    }

    // TODO remove rings
    std::copy(container.begin(), container.end(),
              rings.bound_vortices().begin());
    OutputSnapshot(i, t, rings);

    // 放出する渦を求める
    // TODO 変形したときに位置が合わない
    // TODO remove rings
    auto edge_first = container.edge_begin();
    auto edge_last = container.edge_end();
    std::vector<UVLM::VortexRing> shed(std::distance(edge_first, edge_last));
    UVLM::ShedAtTrailingEdge(edge_first, edge_last,
                             std::begin(shed), rings, Vinfty, t, dt);
    // UVLM::ShedAtTrailingEdge(edge_first, edge_last, shed.begin(),
    //                          vortices->cbegin(), vortices->cend(), Vinfty, t,
    //                          dt);
    std::cerr << "@@" << shed.size() << "\n";

    // TODO remove rings
    UVLM::AdvectWake(&rings, Vinfty, dt);
    // UVLM::AdvectWake(vortices->begin() + wake_offset, vortices->end(),
    //                  vortices->cbegin(), vortices->cend(), Vinfty, dt);
    // UVLM::AdvectWake(rings.wake_vortices().begin(), rings.wake_vortices().end(),
    //                  vortices->cbegin(), vortices->cend(), Vinfty, dt);
    std::vector<UVLM::VortexRing> advected_wake = rings.wake_vortices();

    std::cerr << container.size() << "\n";
    std::cerr << vortices->size() - wake_offset << " vs " << rings.wake_vortices().size() << ">\n";

    // 変形する
    for (auto& vortex : container) {
      for (std::size_t i=0; i<vortex.nodes().size(); i++) {
        morphing.Perfome(&vortex.nodes()[i], vortex.nodes0()[i], t, dt);
      }
    }

    // 変形後のedgeの位置とshedを合わせる
    UVLM::AttachShedVorticesToEdge(container.edge_begin(), container.edge_end(),
                                   shed.begin());

    // 放出した渦の追加
    vortices->insert(vortices->end(), shed.cbegin(), shed.cend());

    // 放出した渦を追加する 
    // rings.AppendWake(shed.cbegin(), shed.cend());

    // TODO remove rings
    std::copy(advected_wake.begin(), advected_wake.end(),
              vortices->begin() + container.size());
    rings.wake_vortices().resize(vortices->size() - wake_offset);
    std::copy(vortices->begin() + wake_offset, vortices->end(),
              rings.wake_vortices().begin());
    for (std::size_t i = 0; i < container.size(); i++) {
      std::cerr << vortices->at(i).gamma() << " vs " << container[i].gamma()
                << " vs " << rings.bound_vortices()[i].gamma() << std::endl;
    }
  }
}

int main(int argc, char* argv[]) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  SimulationBody();
  return 0;
}

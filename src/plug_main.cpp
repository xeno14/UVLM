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

#include <gflags/gflags.h>
#include <iostream>
#include <fstream>
#include <cstdlib>

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

void Output(std::ofstream& ofs, const UVLM::UVLMVortexRing& rings) {
  ofs << std::scientific;
  for (const auto& v : rings.bound_vortices()) {
    for (const auto& vertex : v.nodes()) {
      ofs << vertex.x() << "\t" << vertex.y() << "\t" << vertex.z() << "\t"
          << v.gamma() << std::endl;
    }
  }
  for (const auto& v : rings.wake_vortices()) {
    for (const auto& vertex : v.nodes()) {
      ofs << vertex.x() << "\t" << vertex.y() << "\t" << vertex.z() << "\t"
          << v.gamma() << std::endl;
    }
  }
  ofs << std::endl << std::endl;
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
  UVLM::UVLMVortexRing rings;
  UVLM::Morphing morphing;
  Eigen::Vector3d Vinfty(2, 0, 0.1);

  const Eigen::Vector3d origin(0, 1, 0);
  morphing.set_origin(origin);
  rings.set_origin(origin);

  InitWing(&rings);

  // morphing.set_plug([](double t) { return 0.2 * sin(5*t); });
  morphing.set_flap([](double t) { return M_PI/6 * sin(4*t); });

  // Output(ofs, rings);
  const double dt = FLAGS_dt;

  // main loop
  for(std::size_t i=0; i<100; i++) {
    std::cerr << i << std::endl;
    const double t = i * dt;

    // 連立方程式を解いて翼の上の循環を求める
    auto gamma = UVLM::SolveLinearProblem(
        rings.bound_vortices().cbegin(), rings.bound_vortices().cend(),
        rings.wake_vortices().cbegin(), rings.wake_vortices().cend(), Vinfty,
        morphing, t);
    // std::cerr << gamma << std::endl;
    for (std::size_t i=0; i < rings.bound_vortices().size(); i++) {
      rings.bound_vortices()[i].set_gamma(gamma(i));
    }
    // Output(ofs, rings);
    OutputSnapshot(i, t, rings);

    // 放出する渦を求める
    // TODO 変形したときに位置が合わない
    auto trailing_edge = rings.TrailingEdgeIterators();
    std::vector<UVLM::VortexRing> shed(trailing_edge.second - trailing_edge.first);
    UVLM::ShedAtTrailingEdge(trailing_edge.first, trailing_edge.second,
                             std::begin(shed), rings, morphing, Vinfty, t, dt);

    // 後流の移流
    UVLM::AdvectWake(&rings, Vinfty, dt);

    // 放出した渦を追加する 
    rings.AppendWake(std::begin(shed), std::end(shed));

    // 変形する
    // TODO 関数にする
    for (auto& w : rings.bound_vortices()) {
      for (int i=0; i<4; i++) {
        morphing.Perfome(&w.nodes()[i], w.nodes0()[i], t, dt);
      }
    }
  }
}

int main(int argc, char* argv[]) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  SimulationBody();
  return 0;
}

/**
 * @file plug_main.cpp
 * @brief Add description here
 */

#include "uvlm_vortex_ring.h"
#include "morphing.h"
#include "linear.h"
#include "shed.h"

#include <gflags/gflags.h>
#include <iostream>
#include <fstream>
#include <cstdlib>

DEFINE_string(output, "", "output path");
//DEFINE_string(format, "tsv", "file format");
DEFINE_string(wing_tsv, "", "tsv file");
DEFINE_int32(cols, 0, "number of columns of vertices");
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
  auto pos = ReadWingTsv(FLAGS_wing_tsv);
  rings->InitWing(pos, FLAGS_cols);
}

void Output(std::ofstream& ofs, const UVLM::UVLMVortexRing& rings) {
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

void SimulationBody() {
  UVLM::UVLMVortexRing rings;
  UVLM::Morphing morphing;  // do nothing
  Eigen::Vector3d Vinfty(1, 0, 0.1);

  InitWing(&rings);

  std::ofstream ofs(FLAGS_output);
  if (!ofs) {
    std::cerr << "output open error" << std::endl; 
    std::exit(EXIT_FAILURE);
  }


  Output(ofs, rings);
  const double dt = FLAGS_dt;

  // main loop
  for(int i=0; i<100; i++) {
    std::cerr << i << std::endl;

    const double t = i * dt;
    auto gamma = UVLM::SolveLinearProblem(rings, Vinfty, morphing, t);
    std::cerr << gamma << std::endl;
    for (std::size_t i=0; i < rings.bound_vortices().size(); i++) {
      rings.bound_vortices()[i].set_gamma(gamma(i));
    }
    Output(ofs, rings);
    return;
    auto trailing_edge = rings.TrailingEdgeIterators();
    std::vector<UVLM::VortexRing> shed(trailing_edge.second - trailing_edge.first);
    UVLM::ShedAtTrailingEdge(trailing_edge.first, trailing_edge.second,
                             std::begin(shed), rings, dt);

    rings.AppendWake(std::begin(shed), std::end(shed));
    Output(ofs, rings);
  }
}

int main(int argc, char* argv[]) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  SimulationBody();
  return 0;
}

/**
 * @file result_to_paraview.cpp
 * @brief Add description here
 */


#include "util.h"
#include "../proto/uvlm.pb.h"

#include <gflags/gflags.h>
#include <cstdio>
#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>


DEFINE_string(input, "", "input files with format. e.g. \%08d.dat");
DEFINE_string(output, "stdout", "output destination [stdout, path to file]");
DEFINE_string(snapshot_version, "2", "input file is snapshot1 or snapshot2");

using namespace UVLM;


template <class InputIterator, class F>
void ForEachVortexRingImpl(InputIterator first, InputIterator last, F&& op) {
  while (first != last) {
    op(*first);
    ++first;
  }
}

template <class F>
void ForEachVortexRing(const proto::FlyingWing& wing, F&& op) {
  ForEachVortexRingImpl(wing.bound_vortices().begin(),
                        wing.bound_vortices().end(), op);
  ForEachVortexRingImpl(wing.wake_vortices().begin(),
                        wing.wake_vortices().end(), op);
}

void PlotFlyingWing(FILE* fp, const proto::FlyingWing& wing) {
  ForEachVortexRing(wing, [fp](const auto& ring) {
    for (const auto& node : ring.nodes()) {
      fprintf(fp, "%e %e %e\n", node.x(), node.y(), node.z()); 
    }
  });
}

void PlotGamma(FILE* fp, const proto::FlyingWing& wing) {
  ForEachVortexRing(wing, [fp](const auto& ring) {
    fprintf(fp, "%e\n", ring.gamma());
  });
}

template <class Range>
std::size_t CountTotalSize(const Range& flying_wings) {
  std::size_t res = 0;
  for (const auto& wing : flying_wings) {
    res += wing.bound_vortices().size();
    res += wing.wake_vortices().size();
  }
  return res;
}


template <class Range>
void FlyingWingsToVtk(FILE* fp, const Range& flying_wings) {
  // vortex ringの総数
  const std::size_t total_size = CountTotalSize(flying_wings); 

  // POINTS
  fprintf(fp, "POINTS %lu float\n", total_size * 4);
  for (const auto& wing : flying_wings) {
    PlotFlyingWing(fp, wing);
  }
  fprintf(fp, "\n");

  // CELLS
  fprintf(fp, "CELLS %lu %lu\n", total_size, total_size * 5);
  for (std::size_t i = 0; i < total_size; ++i) {
    fprintf(fp, "4 %lu %lu %lu %lu\n", i*4, i*4+1, i*4+2, i*4+3);
  }
  fprintf(fp, "\n");

  // CELL_TYPES
  fprintf(fp, "CELL_TYPES %lu\n", total_size);
  for (std::size_t i = 0; i < total_size; ++i) {
    fprintf(fp, "9\n");
  }

  // CELL_DATA
  fprintf(fp, "CELL_DATA %lu\n", total_size);
  fprintf(fp, "SCALARS cell_scalars float\n");
  fprintf(fp, "LOOKUP_TABLE default\n");
  for (const auto& wing : flying_wings) {
    PlotGamma(fp, wing);
  }
  fprintf(fp, "\n");
}

void TransformSnapshot2(FILE* fp, const proto::Snapshot2& snapshot) {
  const std::size_t total_size = snapshot.vortices().size();

  // POINTS
  fprintf(fp, "POINTS %lu float\n", total_size * 4);
  for (const auto& vortex : snapshot.vortices()) {
    for (const auto& node : vortex.nodes()) {
      fprintf(fp, "%e %e %e\n", node.x(), node.y(), node.z()); 
    }
  }

  // CELLS
  fprintf(fp, "CELLS %lu %lu\n", total_size, total_size * 5);
  for (std::size_t i = 0; i < total_size; ++i) {
    fprintf(fp, "4 %lu %lu %lu %lu\n", i*4, i*4+1, i*4+2, i*4+3);
  }
  fprintf(fp, "\n");

  // CELL_TYPES
  fprintf(fp, "CELL_TYPES %lu\n", total_size);
  for (std::size_t i = 0; i < total_size; ++i) {
    fprintf(fp, "9\n");
  }

  // CELL_DATA
  fprintf(fp, "CELL_DATA %lu\n", total_size);
  fprintf(fp, "SCALARS cell_scalars float\n");
  fprintf(fp, "LOOKUP_TABLE default\n");
  for (const auto& vortex : snapshot.vortices()) {
    fprintf(fp, "%e\n", vortex.gamma());
  }
  fprintf(fp, "\n");
}

int main(int argc, char* argv[]) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  
  FILE* fp = nullptr;
  if (FLAGS_output == "stdout") {
    fp = stdout;
  } else {
    fp = fopen(FLAGS_output.c_str(), "w");
  }
  CHECK_OPEN(fp);
  std::cerr << "Output: " << FLAGS_output << std::endl;

  // header
  fprintf(fp, "# vtk DataFile Version 2.0\n");
  fprintf(fp, "quad\n");
  fprintf(fp, "ASCII\n");
  fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

  std::ifstream ifs(FLAGS_input, std::ios::binary);
  CHECK_OPEN(ifs);
  std::cerr << "Input : " << FLAGS_input << std::endl;

  if (FLAGS_snapshot_version == "1") {
    proto::Snapshot snapshot;
    snapshot.ParseFromIstream(&ifs);
    // body
    FlyingWingsToVtk(fp, snapshot.flying_wings());
  } else if (FLAGS_snapshot_version == "2") {
    proto::Snapshot2 snapshot;
    snapshot.ParseFromIstream(&ifs);
    TransformSnapshot2(fp, snapshot);
  } else {
    std::cerr << "Invalid version: " << FLAGS_snapshot_version << std::endl;
  }

  return 0;
}

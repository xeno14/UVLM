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

using namespace UVLM;

template <class InputIterator>
void PlotVortexRings(FILE* fp, InputIterator first, InputIterator last) {
  while(first != last) {
    for (const auto& node : first->nodes()) {
      fprintf(fp, "%e %e %e\n", node.x(), node.y(), node.z()); 
    }
    ++first;
  }
}

void PlotFlyingWing(FILE* fp, const proto::FlyingWing& wing) {
  PlotVortexRings(fp, wing.bound_vortices().begin(),
                  wing.bound_vortices().end());
  PlotVortexRings(fp, wing.wake_vortices().begin(), wing.wake_vortices().end());
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

  std::ifstream ifs(FLAGS_input, std::ios::binary);
  CHECK_OPEN(ifs);
  std::cerr << "Input : " << FLAGS_input << std::endl;
  proto::Snapshot snapshot;
  snapshot.ParseFromIstream(&ifs);

  fprintf(fp, "# vtk DataFile Version 2.0\n");
  fprintf(fp, "quad\n");
  fprintf(fp, "ASCII\n");
  fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
  FlyingWingsToVtk(fp, snapshot.flying_wings());

  return 0;
}


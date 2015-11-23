/**
 * @file result_to_paraview.cpp
 * @brief Convert snaphot to vtk file
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

namespace {

template <class InputIterator, class F>
void ForEachVortexRingImpl(InputIterator first, InputIterator last, F&& op) {
  while (first != last) {
    op(*first);
    ++first;
  }
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

void WriteHeader(FILE* fp) {
  fprintf(fp, "# vtk DataFile Version 2.0\n");
  fprintf(fp, "quad\n");
  fprintf(fp, "ASCII\n");
  fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
}

}  // anonymous namespace


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

  WriteHeader(fp);

  std::ifstream ifs(FLAGS_input, std::ios::binary);
  CHECK_OPEN(ifs);
  std::cerr << "Input : " << FLAGS_input << std::endl;

  proto::Snapshot2 snapshot;
  snapshot.ParseFromIstream(&ifs);
  TransformSnapshot2(fp, snapshot);

  return 0;
}

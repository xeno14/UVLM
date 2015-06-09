/**
 * @file result_viewer.cpp
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
    if (first->nodes().size() > 0) {
      const auto& node = *(first->nodes().begin());
      fprintf(fp, "%e %e %e\n", node.x(), node.y(), node.z()); 
    }
    ++first;
    fprintf(fp, "\n\n");
  }
}

void PlotFlyingWing(FILE* fp, const proto::FlyingWing& wing) {
  PlotVortexRings(fp, wing.bound_vortices().begin(),
                  wing.bound_vortices().end());
  PlotVortexRings(fp, wing.wake_vortices().begin(), wing.wake_vortices().end());
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

  fprintf(fp, "splot \"-\" usi 1:2:3 w l title \"t=%e\"\n", snapshot.t());
  for (const auto& wing : snapshot.flying_wings()) {
    PlotFlyingWing(fp, wing);
  }
  fprintf(fp, "end\n");
  fprintf(fp, "\n");

  return 0;
}

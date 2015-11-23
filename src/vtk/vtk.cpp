/**
 * @file vtk.cpp
 * @brief Add description here
 */

#include "vtk.h"

namespace UVLM {
namespace vtk {

void WriteSnapshot2(FILE* fp, const proto::Snapshot2& snapshot) {
  fprintf(fp, "# vtk DataFile Version 2.0\n");
  fprintf(fp, "quad\n");
  fprintf(fp, "ASCII\n");
  fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

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
 
}  // namespace vtk
}  // namespace UVLM

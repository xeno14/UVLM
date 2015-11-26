/**
 * @file velocity.h
 * @brief Add description here
 */
#pragma once

#include "multiple_sheet/multiple_sheet.h"
#include "vortex.h"

using multiple_sheet::MultipleSheet;

namespace UVLM {

inline Eigen::Vector3d InducedVelocity(
    const Eigen::Vector3d& x, const MultipleSheet<Eigen::Vector3d>& pos,
    const MultipleSheet<double>& gamma) {
  Eigen::Vector3d res = Eigen::Vector3d::Zero();
  for (auto index : gamma.list_index()) {
    std::size_t n, i, j;
    std::tie(std::ignore, n, i, j) = index;
    res += UVLM::VORING(x, pos, gamma, n, i, j);
  }
  return res;
}

inline Eigen::Vector3d ParallelInducedVelocity(
    const Eigen::Vector3d& x, const MultipleSheet<Eigen::Vector3d>& pos,
    const MultipleSheet<double>& gamma) {
  Eigen::Vector3d res = Eigen::Vector3d::Zero();
  double vx = 0, vy = 0, vz = 0;
  const auto indices = gamma.list_index();
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : vx, vy, vz)
#endif
  for (std::size_t K = 0; K < indices.size(); ++K) {
    Eigen::Vector3d v;
    std::size_t n, i, j;
    std::tie(std::ignore, n, i, j) = indices[K];
    v = UVLM::VORING(x, pos, gamma, n, i, j);
    vx += v.x();
    vy += v.y();
    vz += v.z();
  }
  return res;
}

}  // namespace UVLM

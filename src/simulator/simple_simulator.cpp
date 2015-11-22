/**
 * @file simple_simulator.cpp
 * @brief Add description here
 */

#include "simple_simulator.h"

namespace UVLM {
namespace simulator {

std::vector<Eigen::Vector3d> CollocationPoints(
    const MultipleSheet<Eigen::Vector3d>& pos) {
  std::vector<Eigen::Vector3d> res;
  for (std::size_t n = 0; n < pos.num(); n++) {
    for (std::size_t i = 0; i < pos.rows() - 1; i++) {
      for (std::size_t j = 0; j < pos.cols() - 1; j++) {
        res.push_back((pos.at(n, i, j) + pos.at(n, i + 1, j) +
                       pos.at(n, i + 1, j + 1) + pos.at(n, i, j + 1)) /
                      4);
      }
    }
  }
  return res;
}

std::vector<Eigen::Vector3d> Normals(
    const MultipleSheet<Eigen::Vector3d>& pos) {
  std::vector<Eigen::Vector3d> res;
  for (std::size_t n = 0; n < pos.num(); n++) {
    for (std::size_t i = 0; i < pos.rows() - 1; ++i) {
      for (std::size_t j = 0; j < pos.cols() - 1; ++j) {
        Eigen::Vector3d nrml =
            (pos.at(n, i + 1, j + 1) - pos.at(n, i, j))
                .cross(pos.at(n, i, j + 1) - pos.at(n, i + 1, j));
        nrml.normalize();
        res.emplace_back(nrml);
      }
    }
  }
  return res;
}

}  // namespace simulator
}  // namespace UVLM

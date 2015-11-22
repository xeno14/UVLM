/**
 * @file simple_simulator.cpp
 * @brief Add description here
 */

#include "simple_simulator.h"

#include <glog/logging.h>

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

void SimpleSimulator::MainLoop(const std::size_t step, const double dt) {
  const double t = step * dt;

  std::copy(wing_gamma.begin(), wing_gamma.end(), wing_gamma_prev.begin());
  
  LOG(INFO) << "Morphing";
  CHECK(wing_pos.size() == wing_pos_init.size());
  for (std::size_t n = 0; n < wing_pos.num(); n++) {
    std::transform(wing_pos_init.sheet_begin(n), wing_pos_init.sheet_end(n),
                   wing_pos.sheet_begin(n),
                   [t, n, this](const Eigen::Vector3d& x0) {
                     return this->morphings_[n].Perfome(x0, t);
                   });
  }
}

}  // namespace simulator
}  // namespace UVLM

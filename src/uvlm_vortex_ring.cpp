/**
 * @file uvlm_vortex_ring.cpp
 * @brief Add description here
 */

#include <iostream>
#include "uvlm_vortex_ring.h"


namespace UVLM {
namespace internal {

void InducedVelocityByVortices(Eigen::Vector3d* const result,
                               const Eigen::Vector3d& pos,
                               const std::vector<VortexRing>& vortices) {
  *result = Eigen::Vector3d::Zero();
  Eigen::Vector3d tmp;
  for (const auto& vortex : vortices) {
    vortex.BiotSavartLaw(&tmp, pos);
    *result += tmp;
  }
}

}  // namespace internal

void UVLMVortexRing::InitWing(const std::vector<Eigen::Vector3d>& pos,
                              std::size_t cols) {
  const std::size_t cols_pos = cols + 1;
  cols_ = cols;
  rows_ = pos.size() / cols_pos - 1;

  for (std::size_t i=0; i<rows_; i++) {  // loop for x
    // y >= 0
    for (std::size_t j=0; j<cols_; j++) {  // loop for y
      std::size_t indices[4] = {
          j + i * cols_pos,            // 0
          j + (i + 1) * cols_pos,      // 1
          j + 1 + (i + 1) * cols_pos,  // 2
          j + 1 + i * cols_pos,        // 3
      };
      VortexRing vortex;
      for (auto idx : indices) vortex.nodes().push_back(pos[idx]);
      bound_vortices_.push_back(vortex);
    }
    // y <= 0
    for (std::size_t j = 0; j < cols_; j++) {  // loop for y
      bound_vortices_.emplace_back();
      bound_vortices_.rbegin()->nodes().resize(4);
    }
  }
  PlaneSymmetry();
  for (auto& vortex : bound_vortices_) vortex.SaveReferenceNode();
}

void UVLMVortexRing::InducedVelocity(Eigen::Vector3d* const result,
                                     const Eigen::Vector3d& pos) const {
  Eigen::Vector3d v, w;
  InducedVelocityByBound(&v, pos);
  InducedVelocityByWake(&w, pos);
  *result = v + w;
}

void UVLMVortexRing::PlaneSymmetry() {
  for (std::size_t i=0; i<rows(); i++) {
    for (std::size_t j=0; j<cols(); j++) {
      const std::size_t idx = j + i * 2 * cols();
      auto& to = bound_vortices_[idx + cols()].nodes();
      const auto& from = bound_vortices_[idx].nodes();
      internal::CopySymmetry(&to[0], from[3]);
      internal::CopySymmetry(&to[1], from[2]);
      internal::CopySymmetry(&to[2], from[1]);
      internal::CopySymmetry(&to[3], from[0]);
    }
  }
}

}  // namespace UVLM

/**
 * @file uvlm_vortex_ring.cpp
 * @brief Add description here
 */

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
  rows_ = pos.size() / cols - 1;

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
      vortex.SaveReferenceNode();
      bound_vortices_.push_back(vortex);
    }
    // y <= 0
    for (std::size_t j = 0; j < cols_; j++) {  // loop for y
      VortexRing vortex;
      const std::size_t idx = j + i * 2 * cols_;  // 対応する渦輪のindex

      // y → -yの変換
      auto f = [](const auto& p) {
        return Eigen::Vector3d(p.x(), -p.y(), p.z());
      };

      // 順番が時計回りになるように入れる
      vortex.nodes().push_back(f(bound_vortices_[idx].nodes()[3]));
      vortex.nodes().push_back(f(bound_vortices_[idx].nodes()[2]));
      vortex.nodes().push_back(f(bound_vortices_[idx].nodes()[1]));
      vortex.nodes().push_back(f(bound_vortices_[idx].nodes()[0]));
      vortex.SaveReferenceNode();
      bound_vortices_.push_back(vortex);
    }
  }
}

void UVLMVortexRing::InducedVelocity(Eigen::Vector3d* const result,
                                     const Eigen::Vector3d& pos) const {
  Eigen::Vector3d v, w;
  InducedVelocityByBound(&v, pos);
  InducedVelocityByWake(&w, pos);
  *result = v + w;
}

}  // namespace UVLM

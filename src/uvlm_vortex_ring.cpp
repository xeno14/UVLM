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

void UVLMVortexRing::InducedVelocity(Eigen::Vector3d* const result,
                                     const Eigen::Vector3d& pos) const {
  Eigen::Vector3d v, w;
  InducedVelocityByBound(&v, pos);
  InducedVelocityByWake(&w, pos);
  *result = v + w;
}

}  // namespace UVLM

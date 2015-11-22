
/**
 * @file calc_impulse.h
 * @brief Add description here
 */
#pragma once

#include "util.h"
#include "uvlm_vortex_ring.h"

namespace UVLM {
namespace calc_impulse {

inline Eigen::Vector3d CalcLambVectorOnPanel(const UVLM::VortexRing& v,
    const std::vector<Eigen::Vector3d>& v_nodes,
    const Eigen::Vector3d& freestream,
    std::vector<VortexRing> vortices) {
  Eigen::Vector3d res = Eigen::Vector3d::Zero();
  auto um = v_nodes.begin();
  v.ForEachSegment([&](const auto& start, const auto& end) {
    const Eigen::Vector3d dl = end - start;
    const Eigen::Vector3d pos = (start + end) / 2;
    Eigen::Vector3d u;
    UVLM::internal::InducedVelocityByVortices(&u, pos, vortices);
    const Eigen::Vector3d ue = *um + u - freestream;
    res += ue.cross(dl) * v.gamma();
    ++um;
  });
  return res;
}

}  // namespace calc_impulse
}  // namespace UVLM

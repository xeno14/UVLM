/**
 * @file shed.cpp
 * @brief Add description here
 */

#include "shed.h"

#include <iostream>

namespace UVLM {
namespace internal {

void ShedSingleAtTrailingEdge(VortexRing* result, const VortexRing& target,
                              const UVLMVortexRing& rings,
                              const Eigen::Vector3d& Vinfty, const double t,
                              const double dt) {
  *result = target;

  Eigen::Vector3d velocity;
  for (auto& node : result->nodes()) {
    rings.InducedVelocity(&velocity, node);
    velocity += Vinfty;
    Advect(&node, velocity, dt);
  }
}

void AdvectWakeImpl(std::vector<VortexRing>* result,
                    const UVLMVortexRing& rings, const Eigen::Vector3d& Vinfty,
                    const double dt) {
  for (auto& wake : *result) {
    Eigen::Vector3d velocity;
    for (auto& node : wake.nodes()) {
      rings.InducedVelocity(&velocity, node);
      velocity += Vinfty;
      internal::Advect(&node, velocity, dt);
    }
  }
}

}  // namespace internal

void AdvectWake(UVLMVortexRing* rings, const Eigen::Vector3d& Vinfty,
                const double dt) {
  std::vector<VortexRing> new_wake;
  new_wake = rings->wake_vortices();
  internal::AdvectWakeImpl(&new_wake, *rings, Vinfty, dt);
  rings->wake_vortices().swap(new_wake);
}

}  // namespace UVLM

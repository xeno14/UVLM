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
  const std::size_t sz = result->size();
  std::size_t i = 0;
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (i=0; i < sz; i++) {
    auto& wake = result->at(i);
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

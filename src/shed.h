
/**
 * @file shed.h
 * @brief Add description here
 */
#pragma once

#include "vortex.h"
#include "uvlm_vortex_ring.h"

namespace UVLM {
namespace internal {

void Advect(Eigen::Vector3d* target, const Eigen::Vector3d& vel,
            const double dt);

void ShedSingleAtTrailingEdge(VortexRing* result, const VortexRing& target,
                              const UVLMVortexRing& rings,
                              const Eigen::Vector3d& Vinfty, const double dt);

}  // namespace internal

template <class InputIterator, class OutputIterator>
void ShedAtTrailingEdge(InputIterator first, InputIterator last,
                        OutputIterator result, const UVLMVortexRing& rings,
                        const Eigen::Vector3d& Vinfty, const double dt) {
  while (first != last) {
    internal::ShedSingleAtTrailingEdge(&(*result), *first, rings, Vinfty, dt);
    ++first; ++result;
  }
}

void AdvectWake(UVLMVortexRing* rings, const Eigen::Vector3d& Vinfty,
                const double dt);

}  // namespace UVLM

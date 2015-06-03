
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
                              const UVLMVortexRing& rings, const double dt);

}  // namespace internal

template <class InputIterator, class OutputIterator>
void ShedAtTrailingEdge(InputIterator first, InputIterator last,
                        OutputIterator result, const UVLMVortexRing& rings,
                        const double dt) {
  while (first != last) {
    internal::ShedSingleAtTrailingEdge(&(*result), *first, rings, dt);
    ++first; ++result;
  }
}

}  // namespace UVLM

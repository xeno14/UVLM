
/**
 * @file shed.h
 * @brief Add description here
 */
#pragma once

#include "vortex.h"

namespace UVLM {
namespace internal {

void ShedSingleVortex(VortexRing* result, const VortexRing& target,
                      const std::vector<VortexRing>& vortices,
                      const std::vector<VortexRing>& wake, double dt);

}  // namespace internal

void Advect(Eigen::Vector3d* target, const Eigen::Vector3d& vel,
            const double dt);

void InducedVelocity(Eigen::Vector3d* const result, const Eigen::Vector3d& pos,
              const std::vector<VortexRing>& vortices);

void InducedVelocity(Eigen::Vector3d* const result, const Eigen::Vector3d& pos,
              const std::vector<VortexRing>& vortices,
              const std::vector<VortexRing>& wake);

}  // namespace UVLM

/**
 * @file shed.cpp
 * @brief Add description here
 */

#include "shed.h"


namespace UVLM {
namespace internal {

void ShedSingleVortex(VortexRing* result, const VortexRing& target,
                      const std::vector<VortexRing>& vortices,
                      const std::vector<VortexRing>& wake, double dt) {
  //
  // 3--2=3'---2'
  // |   |     |
  // 0--1=0'---1'
  result->nodes()[0] = target.nodes()[1];
  result->nodes()[3] = target.nodes()[2];
  result->nodes()[1] = target.nodes()[1];
  result->nodes()[2] = target.nodes()[2];

  Eigen::Vector3d v1, v2;
  InducedVelocity(&v1, target.nodes()[1], vortices, wake);
  InducedVelocity(&v2, target.nodes()[2], vortices, wake);
  Advect(&result->nodes()[1], v1, dt);
  Advect(&result->nodes()[2], v2, dt);
}

}  // namespace internal

void Advect(Eigen::Vector3d* target, const Eigen::Vector3d& vel,
            const double dt) {
  *target += vel * dt;  
}

void InducedVelocity(Eigen::Vector3d* const result, const Eigen::Vector3d& pos,
              const std::vector<VortexRing>& vortices) {
  *result << 0, 0, 0;
  Vector3d tmp;
  for (const auto& vortex : vortices) {
    vortex.BiotSavartLaw(&tmp, pos);
    *result += tmp;
  }
}

void InducedVelocity(Eigen::Vector3d* const result, const Eigen::Vector3d& pos,
              const std::vector<VortexRing>& vortices,
              const std::vector<VortexRing>& wake) {
  Eigen::Vector3d v, w;
  InducedVelocity(&v, pos, vortices);
  InducedVelocity(&w, pos, wake);
  *result = v + w;
}

}  // namespace UVLM

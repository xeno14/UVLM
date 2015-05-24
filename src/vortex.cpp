/**
 * @file vortex.cpp
 * @brief Add description here
 */

#include "vortex.h"

namespace UVLM {

double DistanceLineAndPoint(const Vector3d& start,
                            const Vector3d& end,
                            const Vector3d& pos) {
  Vector3d direction = end - start;
  direction.normalize();

  Vector3d diff = pos - start;
  return (diff - direction.dot(diff) * direction).norm();
}

void BiotSavartLaw(Vector3d* result,
                   const Vector3d& start, const Vector3d& end,
                   const Vector3d& pos) {
  double h = DistanceLineAndPoint(start, end, pos);

  Vector3d direction = end - start;
  direction.normalize();

  Vector3d diff1 = pos - start; diff1.normalize();
  Vector3d diff2 = pos - end;   diff2.normalize();
  double cos1 = diff1.dot(direction);
  double cos2 = diff2.dot(direction);

  *result = direction.cross(diff1);
  result->normalize();
  *result *= (cos1 - cos2) / (4. * M_PI * h);
}

void VortexFilament::BiotSavartLaw(Vector3d* result,
                                const Vector3d& pos) const {
    ::UVLM::BiotSavartLaw(result, start_, end_, pos);
    *result *= gamma_;
}

void VortexRing::BiotSavartLaw(Vector3d* result, const Vector3d& pos) const {
  *result << 0, 0, 0;
  Vector3d tmp;
  for (const auto& filament : filaments_) {
    filament.BiotSavartLaw(&tmp, pos); *result += tmp;
  }
}

VortexRing& VortexRing::PushNode(const Vector3d& pos) {
  nodes_.push_back(pos);
  return *this;
}

void VortexRing::AssembleRing() {
  if (nodes_.size() < 3) return; // exception?

  nodes_.push_back(nodes_[0]);  // for circulation
  for (std::size_t i=0; i<nodes_.size() -1; i++) {
    filaments_.emplace_back(gamma_, nodes_[i], nodes_[i+1]);
  }
  nodes_.clear();
}

}  // namespace UVLM

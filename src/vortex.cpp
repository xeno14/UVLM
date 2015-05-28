/**
 * @file vortex.cpp
 * @brief Add description here
 */

#include "vortex.h"

#include <algorithm>

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

void VortexRing::BiotSavartLaw(Vector3d* result, const Vector3d& pos,
                               double gamma) const {
  *result << 0, 0, 0;
  Vector3d tmp;
  for (std::size_t i = 0; i < nodes_.size(); i++) {
    const Vector3d& start = nodes_[i];
    const Vector3d& end = nodes_[(i + 1) % nodes_.size()];
    ::UVLM::BiotSavartLaw(&tmp, start, end, pos);
    *result += tmp;
  }
  *result *= gamma;
}

VortexRing& VortexRing::PushNode(const Vector3d& pos) {
  nodes_.push_back(pos);
  return *this;
}

VortexRing& VortexRing::PushNode(double x, double y, double z) {
  nodes_.emplace_back(x, y, z);
  return *this;
}

void VortexRing::SaveReferenceNode() {
  nodes0_.resize(nodes_.size());
  std::copy(std::begin(nodes_), std::end(nodes_), std::begin(nodes0_));
}

Vector3d VortexRing::Normal() const {
  Vector3d diff1 = nodes_[2] - nodes_[0];
  Vector3d diff2 = nodes_[3] - nodes_[1];
  return diff1.cross(diff2);
}

Vector3d VortexRing::Centroid() const {
  // TODO ちゃんと交点を計算する
  // Vector3d res;
  // for (const auto& node : nodes_) {
  //   res += node;
  // }
  // res /= nodes_.size();
  // return res;
  return (nodes_[0] + nodes_[1] + nodes_[2] + nodes_[3]) / 4;
}

}  // namespace UVLM

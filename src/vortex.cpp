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

  // OMAJINAI
  if (std::isnan(result->x()) || std::isnan(result->y()) ||
      std::isnan(result->z()))
    *result = Eigen::Vector3d::Zero();
}

void VortexRing::BiotSavartLaw(Vector3d* result, const Vector3d& pos,
                               double gamma) const {
  *result = ::Eigen::Vector3d::Zero();
  Vector3d tmp;
  this->ForEachSegment(
      [&](const auto& start, const auto& end) {
        ::UVLM::BiotSavartLaw(&tmp, start, end, pos);
        *result += tmp;
      });
  *result *= gamma;
}

void VortexRing::ChordwiseBiotSavartLaw(Vector3d* result, const Vector3d& pos,
                              double gamma) const {
  *result = ::Eigen::Vector3d::Zero();
  Vector3d tmp;
  // 0-1 and 2-3
  static const int CHORDWISE_INDICES[] = {0, 2};
  for (std::size_t i : CHORDWISE_INDICES) {
    const Vector3d& end = nodes_[i];
    const Vector3d& start = nodes_[(i + 1) % nodes_.size()];
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
  Vector3d res = diff1.cross(diff2);
  res.normalize();
  return res;
}

Vector3d VortexRing::Tangent() const {
  Vector3d res = nodes_[2] - nodes_[1] + nodes_[1] - nodes_[3];
  res.normalize();
  return res;
}

Vector3d VortexRing::Centroid() const {
  return (nodes_[0] + nodes_[1] + nodes_[2] + nodes_[3]) / 4;
}

Vector3d VortexRing::ReferenceCentroid() const {
  return (nodes0_[0] + nodes0_[1] + nodes0_[2] + nodes0_[3]) / 4;
}

Eigen::Vector3d VortexRing::TanVecChord() const {
  Eigen::Vector3d res(this->nodes()[1] - this->nodes()[0]);
  res.normalize();
  return res;
}

Eigen::Vector3d VortexRing::TanVecSpan() const {
  Eigen::Vector3d res(this->nodes()[3] - this->nodes()[0]);
  res.normalize();
  return res;
}

double VortexRing::AngleOfAttack(const Eigen::Vector3d& Q) const {
  const double y = Q.dot(Normal());
  const double x = Q.dot(Tangent());
  return atan2(y, x);
}

double VortexRing::CalcC() const {
  return (this->nodes()[0] - this->nodes()[1]).norm();
}

double VortexRing::CalcB() const {
  return (this->nodes()[0] - this->nodes()[3]).norm();
}

}  // namespace UVLM

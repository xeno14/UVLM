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

void BiotSavartLaw(Vector3d* result, const Vector3d& start, const Vector3d& end,
                   const Vector3d& pos) {
  // see Katz and Plotkin p.255
  Eigen::Vector3d r0 = end - start;
  Eigen::Vector3d r1 = pos - start;
  Eigen::Vector3d r2 = pos - end;
  Eigen::Vector3d d = r1.cross(r2);
  double d_len = d.norm();
  double r1_len = r1.norm();
  double r2_len = r2.norm();
  static const double eps = 1e-10;  // cut off length
  if (d_len * d_len < eps || r1_len < eps || r2_len < eps) {
    *result = Eigen::Vector3d::Zero();
    return;
  }
  const double K = 1. / 4. / M_PI / d.squaredNorm() *
                   r0.dot(r1 / r1.norm() - r2 / r2.norm());
  *result = d * K;
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
  // see Ghommem 2011 (5.1)
  Eigen::Vector3d res = (nodes_[1] - nodes_[3]).cross(nodes_[2] - nodes_[0]);
  res.normalize();
  return res;
}

Vector3d VortexRing::Tangent() const {
  // see Ghommem 2011 (5.2)
  Vector3d res = (nodes_[1] - nodes_[3]) + (nodes_[2] - nodes_[0]);
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

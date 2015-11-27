/**
 * @file vortex_kernel.cpp
 * @brief Add description here
 */

#include "vortex_kernel.h"


namespace UVLM {
namespace vortex_kernel {

Eigen::Vector3d CutOffKernel::Induce(const Eigen::Vector3d& x,
                                     const Eigen::Vector3d& x1,
                                     const Eigen::Vector3d& x2,
                                     const double gamma) const {
  const auto r0 = x2 - x1;
  const auto r1 = x - x1;
  const auto r2 = x - x2;
  const auto d = r1.cross(r2);

  if (r1.norm() < cutoff_length_ || r2.norm() < cutoff_length_ ||
      d.squaredNorm() < cutoff_length_) {
    return Eigen::Vector3d::Zero();
  }

  const double coeff = gamma / (4. * M_PI * d.squaredNorm()) *
                       (r0.dot(r1) / r1.norm() - r0.dot(r2) / r2.norm());
  return d * coeff;
}

Eigen::Vector3d RosenheadMooreKernel::Induce(const Eigen::Vector3d& x,
                                             const Eigen::Vector3d& x1,
                                             const Eigen::Vector3d& x2,
                                             const double gamma) const {
  const auto r0 = x2 - x1;
  const auto r1 = x - x1;
  const auto r2 = x - x2;
  const auto e = r1.cross(r2);
  const double d = e.norm() / r0.norm();
  const double c1 = r0.dot(r1)/r0.norm()/r1.norm();
  const double c2 = r0.dot(r2)/r0.norm()/r2.norm();
  const double a = (delta_ / d) * (delta_ / d);

  const double coeff =
      gamma / (4. * M_PI) / (1 + a) / e.squaredNorm() * r0.norm() *
      (c1 / sqrt(1 + a * (1 - c1 * c1)) - c2 / sqrt(1 + a * (1 - c2 * c2)));

  return e * coeff;
}

}  // namespace vortex_kernel
}  // namespace UVLM

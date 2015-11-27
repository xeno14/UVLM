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

}  // namespace vortex_kernel
}  // namespace UVLM

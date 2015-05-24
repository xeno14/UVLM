/**
 * @file vortex.cpp
 * @brief Add description here
 */

#include "vortex.h"

namespace UVLM {

double DistanceLineAndPoint(const Eigen::Vector3d& start,
                            const Eigen::Vector3d& end,
                            const Eigen::Vector3d& pos) {
  Eigen::Vector3d direction = end - start;
  direction.normalize();

  Eigen::Vector3d diff = pos - start;
  return (diff - direction.dot(diff) * direction).norm();
}

void BiotSavartLaw(Eigen::Vector3d* result,
                   const Eigen::Vector3d& start, const Eigen::Vector3d& end,
                   const Eigen::Vector3d& pos) {
  double h = DistanceLineAndPoint(start, end, pos);

  Eigen::Vector3d direction = end - start;
  direction.normalize();

  Eigen::Vector3d diff1 = pos - start; diff1.normalize();
  Eigen::Vector3d diff2 = pos - end;   diff2.normalize();
  double cos1 = diff1.dot(direction);
  double cos2 = diff2.dot(direction);

  *result = direction.cross(diff1);
  result->normalize();
  *result *= (cos1 - cos2) / (4. * M_PI * h);
}

void VortexFilament::BiotSavart(Eigen::Vector3d* result,
                                const Eigen::Vector3d& pos) const {
    BiotSavartLaw(result, start_, end_, pos);
    *result *= gamma_;
}

}  // namespace UVLM

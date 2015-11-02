
/**
 * @file functions.h
 * @brief Add description here
 */
#pragma once

#include "../vortex.h"
#include "../morphing.h"

namespace UVLM {
namespace calc_load {
namespace internal {

inline Eigen::Matrix3d CalcProjectionOperator(const Eigen::Vector3d& Um) {
  Eigen::Vector3d Um_ = Um; // normalized vector
  Um_.normalize();
  return Eigen::Matrix3d::Identity() - Um_ * Um_.transpose();
}

/**
 * @brief Calc velocity contribution from the surface motion
 */
inline Eigen::Vector3d CalcUm(const Morphing& morphing,
                              const Eigen::Vector3d& x0, 
                              const Eigen::Vector3d& freestream,
                              const double t) {
  Eigen::Vector3d Uls;
  morphing.Velocity(&Uls, x0, t);
  return freestream - Uls;
}

}  // namespace internal

struct AerodynamicLoad {
  Eigen::Vector3d F;
  double Pin, Pout;
};

}  // calc_load
}  // UVLM

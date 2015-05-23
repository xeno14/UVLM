/**
 * @file morphing.cpp
 * @brief Add description here
 */

#include "morphing.h"

namespace UVLM {

Morphing::Morphing() {
  alpha_ = 0;
  pluging_ = internal::DefaultFunc;
  flapping_ = internal::DefaultFunc;
}

void Morphing::Perfome(Eigen::Vector3d* x, const Eigen::Vector3d& x0,
                       const double t) const {
  Eigen::Matrix3d T;
  PrepareMatrix(&T, x0, t);

  double gamma_z = pluging_(t);

  *x = Eigen::Vector3d::UnitZ() * gamma_z + T * x0;
}

void Morphing::PrepareMatrix(Eigen::Matrix3d* m, const Eigen::Vector3d& x0,
                             double t) const {
  Eigen::Matrix3d twist;
  twist = Eigen::Matrix3d::Identity();  // TODO

  Eigen::Matrix3d flap;
  double phi = flapping_(t);
  flap << 1, 0, 0,
          0, cos(phi), sin(phi),
          0, -sin(phi), cos(phi);
  Eigen::Matrix3d attack;
  attack << cos(alpha_), 0, sin(alpha_),
            0, 1, 0,
            -sin(alpha_), 0, cos(alpha_);

  *m = attack * flap * twist;
}

}  // namespace UVLM;

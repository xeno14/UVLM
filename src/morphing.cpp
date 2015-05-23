/**
 * @file morphing.cpp
 * @brief Add description here
 */

#include "morphing.h"

namespace UVLM {

Morphing::Morphing() {
  alpha_ = 0;
  plug_ = internal::DefaultFunc;
  flap_ = internal::DefaultFunc;
  twist_ = internal::DefaultFunc2;
}

void Morphing::Perfome(Eigen::Vector3d* x, const Eigen::Vector3d& x0,
                       const double t) const {
  Eigen::Matrix3d T;
  PrepareMatrix(&T, x0, t);

  double gamma_z = plug_(t);

  *x = Eigen::Vector3d::UnitZ() * gamma_z + T * x0;
}

void Morphing::PrepareMatrix(Eigen::Matrix3d* m, const Eigen::Vector3d& x0,
                             double t) const {
  Eigen::Matrix3d twist;
  const double beta = twist_(x0, t);
  twist << cos(beta), 0, sin(beta),
           0, 1, 0,
           -sin(beta), 0, cos(beta);

  Eigen::Matrix3d flap;
  const double phi = flap_(t);
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

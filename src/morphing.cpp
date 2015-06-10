/**
 * @file morphing.cpp
 * @brief Add description here
 */

#include "morphing.h"

namespace UVLM {

Morphing::Morphing() {
  alpha_ = 0;
  Clear();
}

void Morphing::Perfome(Eigen::Vector3d* x, const Eigen::Vector3d& origin,
                       const Eigen::Vector3d& x0, const double t,
                       const double dt) const {
  Eigen::Vector3d x_ref = x0 - origin;
  bool is_negative = x0.y() < 0;

  if (is_negative) x_ref.y() *= -1;

  Eigen::Matrix3d T;
  PrepareMatrix(&T, x_ref, t);

  const double gamma_z = plug_(t);
  const double gamma_b = bend_(x_ref, t);

  *x = Eigen::Vector3d::UnitZ() * gamma_z +
       T * (x_ref + Eigen::Vector3d::UnitZ() * gamma_b);

  if (is_negative) x->y() *= -1;

  *x += origin;
}

void Morphing::Velocity(Eigen::Vector3d* v, const Eigen::Vector3d& origin,
                        const Eigen::Vector3d& x0, const double t,
                        const double dt) const {
  Eigen::Vector3d x1, x2;   // x(t-dt), x(t+dt)
  Perfome(&x1, origin, x0, t - dt, dt);
  Perfome(&x2, origin, x0, t + dt, dt);
  *v = (x2 - x1) / (2*dt);
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

void Morphing::Clear() {
  plug_ =  internal::DefaultFunc;
  flap_ =  internal::DefaultFunc;
  twist_ = internal::DefaultFunc2;
  bend_  = internal::DefaultFunc2;
  // thrust_ = internal::DefaultFunc3;
}

}  // namespace UVLM;

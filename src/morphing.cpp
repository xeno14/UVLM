/**
 * @file morphing.cpp
 * @brief Implemention of morphing.h
 */

#include "morphing.h"

namespace UVLM {

Morphing::Morphing() : origin_(0, 0, 0)  {
  alpha_ = 0;
  Clear();
}

Morphing::Morphing(const Morphing& m)
    : alpha_(m.alpha_),
      plug_(m.plug_),
      flap_(m.flap_),
      twist_(m.twist_),
      bend_(m.bend_),
      origin_(m.origin_) {}

void Morphing::Perfome(Eigen::Vector3d* x, const Eigen::Vector3d& x0,
                       const double t) const {
  Eigen::Vector3d x_ref = x0 - origin_;
  bool is_negative = x_ref.y() < 0;

  if (is_negative) x_ref.y() *= -1;

  Eigen::Matrix3d T;
  PrepareMatrix(&T, x_ref, t);

  const double gamma_z = plug_(t);
  const double gamma_b = bend_(x_ref, t);

  *x = Eigen::Vector3d::UnitZ() * gamma_z +
       T * (x_ref + Eigen::Vector3d::UnitZ() * gamma_b);

  if (x_ref.y() == 0) x->y() = 0;
  if (is_negative) x->y() *= -1;

  *x += origin_;
}

void Morphing::Velocity(Eigen::Vector3d* v,
                        const Eigen::Vector3d& x0, const double t,
                        const double dt) const {
  Eigen::Vector3d x1, x2;   // x(t-dt), x(t+dt)
  Perfome(&x1, x0, t - dt);
  Perfome(&x2, x0, t + dt);
  *v = (x2 - x1) / (2*dt);
}

void Morphing::PrepareMatrix(Eigen::Matrix3d* m, const Eigen::Vector3d& x0,
                             double t) const {
  // see Kats and Plotkin p.371, 423
  // see Ghommem et al. (2012) eq(15)
  Eigen::Matrix3d flap;
  const double phi = flap_(t);
  flap << 1, 0, 0,
          0, cos(phi), sin(phi),
          0, -sin(phi), cos(phi);

  Eigen::Matrix3d attack;
  const double alpha = alpha_;
  attack << cos(alpha), 0, sin(alpha),
            0, 1, 0,
           -sin(alpha), 0, cos(alpha);

  Eigen::Matrix3d twist;
  const double beta = twist_(x0, t);
  twist << cos(beta), 0, sin(beta),
           0,         1, 0,
          -sin(beta), 0, cos(beta);
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

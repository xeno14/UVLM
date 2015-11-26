
/**
 * @file advect.h
 * @brief Add description here
 */
#pragma once

#include "multiple_sheet/multiple_sheet.h"
#include "velocity.h"

using multiple_sheet::MultipleSheet;

namespace UVLM {
namespace advect {

class Advection {
 public:
  Advection() = default;
  virtual ~Advection() = default;
  virtual void Advect(MultipleSheet<Eigen::Vector3d>* next,
                      const MultipleSheet<Eigen::Vector3d>& wing_pos,
                      const MultipleSheet<double>& wing_gamma,
                      const MultipleSheet<Eigen::Vector3d>& wake_pos,
                      const MultipleSheet<double>& wake_gamma,
                      const Eigen::Vector3d& forward_flight,
                      const double dt) const = 0;
 protected:
  Eigen::Vector3d Velocity(const Eigen::Vector3d& x,
                           const MultipleSheet<Eigen::Vector3d>& wing_pos,
                           const MultipleSheet<double>& wing_gamma,
                           const MultipleSheet<Eigen::Vector3d>& wake_pos,
                           const MultipleSheet<double>& wake_gamma,
                           const Eigen::Vector3d& forward_flight) const {
    return -forward_flight + UVLM::InducedVelocity(x, wing_pos, wing_gamma) +
           UVLM::InducedVelocity(x, wake_pos, wake_gamma);
  }
};

class Euler : public Advection {
 public:
  void Advect(MultipleSheet<Eigen::Vector3d>* next,
              const MultipleSheet<Eigen::Vector3d>& wing_pos,
              const MultipleSheet<double>& wing_gamma,
              const MultipleSheet<Eigen::Vector3d>& wake_pos,
              const MultipleSheet<double>& wake_gamma,
              const Eigen::Vector3d& forward_flight,
              const double dt) const override;
 private:
  mutable MultipleSheet<Eigen::Vector3d> vel;
};

class RungeKutta2 : public Advection {
 public:
  void Advect(MultipleSheet<Eigen::Vector3d>* next,
              const MultipleSheet<Eigen::Vector3d>& wing_pos,
              const MultipleSheet<double>& wing_gamma,
              const MultipleSheet<Eigen::Vector3d>& wake_pos,
              const MultipleSheet<double>& wake_gamma,
              const Eigen::Vector3d& forward_flight,
              const double dt) const override;
 private:
  mutable MultipleSheet<Eigen::Vector3d> k1, k2;
  mutable MultipleSheet<Eigen::Vector3d> pos1, pos2;
};

}  // namespace advect
}  // namespace UVLM

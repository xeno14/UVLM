
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

inline Eigen::Vector3d Velocity(const Eigen::Vector3d& x,
                                const MultipleSheet<Eigen::Vector3d>& wing_pos,
                                const MultipleSheet<double>& wing_gamma,
                                const MultipleSheet<Eigen::Vector3d>& wake_pos,
                                const MultipleSheet<double>& wake_gamma,
                                const Eigen::Vector3d& forward_flight) {
  return -forward_flight + UVLM::InducedVelocity(x, wing_pos, wing_gamma) +
         UVLM::InducedVelocity(x, wake_pos, wake_gamma);
}

inline void Velocity(MultipleSheet<Eigen::Vector3d>* result,
                     const MultipleSheet<Eigen::Vector3d>& wing_pos,
                     const MultipleSheet<double>& wing_gamma,
                     const MultipleSheet<Eigen::Vector3d>& wake_pos,
                     const MultipleSheet<double>& wake_gamma,
                     const Eigen::Vector3d& forward_flight) {
  result->resize(wake_pos);
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (std::size_t K = 0; K < wake_pos.size(); K++) {
    (*result)[K] = Velocity(wake_pos[K], wing_pos, wing_gamma, wake_pos,
                            wake_gamma, forward_flight);
  }

  // for trailing edge
  for (std::size_t n = 0; n < result->num(); n++) {
    std::fill(result->iterator_at(n, 0, 0), result->iterator_at(n, 1, 0),
              -forward_flight);
  }
}

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
  RungeKutta2(std::size_t steps = 5) : ignore_steps_(steps) {}
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
  const std::size_t ignore_steps_;
};

class RungeKutta4 : public Advection {
 public:
  RungeKutta4(std::size_t steps = 5) : ignore_steps_(steps) {}
  void Advect(MultipleSheet<Eigen::Vector3d>* next,
              const MultipleSheet<Eigen::Vector3d>& wing_pos,
              const MultipleSheet<double>& wing_gamma,
              const MultipleSheet<Eigen::Vector3d>& wake_pos,
              const MultipleSheet<double>& wake_gamma,
              const Eigen::Vector3d& forward_flight,
              const double dt) const override;

 private:
  mutable MultipleSheet<Eigen::Vector3d> k1, k2, k3, k4;
  mutable MultipleSheet<Eigen::Vector3d> pos1, pos2, pos3, pos4;
  const std::size_t ignore_steps_;
};

class AdamsBashforth2 : public Advection {
 public:
  AdamsBashforth2() {}
  void Advect(MultipleSheet<Eigen::Vector3d>* next,
              const MultipleSheet<Eigen::Vector3d>& wing_pos,
              const MultipleSheet<double>& wing_gamma,
              const MultipleSheet<Eigen::Vector3d>& wake_pos,
              const MultipleSheet<double>& wake_gamma,
              const Eigen::Vector3d& forward_flight,
              const double dt) const override;

 private:
  mutable MultipleSheet<Eigen::Vector3d> v, v_prev;
};

}  // namespace advect
}  // namespace UVLM

/**
 * @file advect.cpp
 * @brief Add description here
 */

#include "advect.h"

namespace UVLM {
namespace advect {

void Euler::Advect(MultipleSheet<Eigen::Vector3d>* next,
                   const vortex_kernel::VortexKernel& kernel,
                   const MultipleSheet<Eigen::Vector3d>& wing_pos,
                   const MultipleSheet<double>& wing_gamma,
                   const MultipleSheet<Eigen::Vector3d>& wake_pos,
                   const MultipleSheet<double>& wake_gamma,
                   const Eigen::Vector3d& forward_flight,
                   const double dt) const {
  Velocity(kernel, &vel, wing_pos, wing_gamma, wake_pos, wake_gamma, forward_flight);
  for (std::size_t K = 0; K < wake_pos.size(); K++) {
    (*next)[K] = wake_pos[K] + vel[K] * dt;
  }
}

void RungeKutta2::Advect(MultipleSheet<Eigen::Vector3d>* next,
                         const vortex_kernel::VortexKernel& kernel,
                         const MultipleSheet<Eigen::Vector3d>& wing_pos,
                         const MultipleSheet<double>& wing_gamma,
                         const MultipleSheet<Eigen::Vector3d>& wake_pos,
                         const MultipleSheet<double>& wake_gamma,
                         const Eigen::Vector3d& forward_flight,
                         const double dt) const {
  k1.resize(wake_pos.num(), wake_pos.rows(), wake_pos.cols());
  k2.resize(wake_pos.num(), wake_pos.rows(), wake_pos.cols());
  pos1.resize(wake_pos.num(), wake_pos.rows(), wake_pos.cols());
  pos2.resize(wake_pos.num(), wake_pos.rows(), wake_pos.cols());

  // 1st step
  for (std::size_t K = 0; K < wake_pos.size(); K++) {
    pos1[K] = wake_pos[K];
  }
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (std::size_t K = 0; K < wake_pos.size(); K++) {
    k1[K] = Velocity(kernel, pos1[K], wing_pos, wing_gamma, pos1, wake_gamma,
                     forward_flight);
  }
  for (std::size_t n = 0; n < k1.num(); n++) {
    // for trailing edge
    std::fill(k1.iterator_at(n, 0, 0), k1.iterator_at(n, 1, 0),
              -forward_flight);
  }

  // 2nd step
  for (std::size_t K = 0; K < wake_pos.size(); K++) {
    pos2[K] = wake_pos[K] + k1[K] * dt;
  }
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (std::size_t K = 0; K < wake_pos.size(); K++) {
    k2[K] = Velocity(kernel, pos2[K], wing_pos, wing_gamma, pos2, wake_gamma,
                     forward_flight);
  }
  for (std::size_t n = 0; n < k2.num(); n++) {
    // for trailing edge
    std::fill(k2.iterator_at(n, 0, 0), k2.iterator_at(n, 1, 0),
              -forward_flight);
  }

  // update
  for (const auto& index : wake_pos.list_index()) {
    std::size_t K, i;
    std::tie(K, std::ignore, i, std::ignore) = index;
    if (i > ignore_steps_) {
      (*next)[K] = wake_pos[K] + (k1[K] + k2[K]) * dt / 2;  // RK2
    } else {
      (*next)[K] = wake_pos[K] + k1[K] * dt;  // Euler
    }
  }
}

void RungeKutta4::Advect(MultipleSheet<Eigen::Vector3d>* next,
                         const vortex_kernel::VortexKernel& kernel,
                         const MultipleSheet<Eigen::Vector3d>& wing_pos,
                         const MultipleSheet<double>& wing_gamma,
                         const MultipleSheet<Eigen::Vector3d>& wake_pos,
                         const MultipleSheet<double>& wake_gamma,
                         const Eigen::Vector3d& forward_flight,
                         const double dt) const {
  k1.resize(wake_pos.num(), wake_pos.rows(), wake_pos.cols());
  k2.resize(wake_pos.num(), wake_pos.rows(), wake_pos.cols());
  k3.resize(wake_pos.num(), wake_pos.rows(), wake_pos.cols());
  k4.resize(wake_pos.num(), wake_pos.rows(), wake_pos.cols());
  pos1.resize(wake_pos.num(), wake_pos.rows(), wake_pos.cols());
  pos2.resize(wake_pos.num(), wake_pos.rows(), wake_pos.cols());
  pos3.resize(wake_pos.num(), wake_pos.rows(), wake_pos.cols());
  pos4.resize(wake_pos.num(), wake_pos.rows(), wake_pos.cols());

  // 1st step
  for (std::size_t K = 0; K < wake_pos.size(); K++) {
    pos1[K] = wake_pos[K];
  }
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (std::size_t K = 0; K < wake_pos.size(); K++) {
    k1[K] = Velocity(kernel, pos1[K], wing_pos, wing_gamma, pos1, wake_gamma,
                     forward_flight);
  }
  for (std::size_t n = 0; n < k1.num(); n++) {
    // for trailing edge
    std::fill(k1.iterator_at(n, 0, 0), k1.iterator_at(n, 1, 0),
              -forward_flight);
  }

  // 2nd step
  for (std::size_t K = 0; K < wake_pos.size(); K++) {
    pos2[K] = wake_pos[K] + k1[K] * dt / 2;
  }
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (std::size_t K = 0; K < wake_pos.size(); K++) {
    k2[K] = Velocity(kernel, pos2[K], wing_pos, wing_gamma, pos2, wake_gamma,
                     forward_flight);
  }
  for (std::size_t n = 0; n < k2.num(); n++) {
    // for trailing edge
    std::fill(k2.iterator_at(n, 0, 0), k2.iterator_at(n, 1, 0),
              -forward_flight);
  }

  // 3rd step
  for (std::size_t K = 0; K < wake_pos.size(); K++) {
    pos3[K] = wake_pos[K] + k2[K] * dt / 2;
  }
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (std::size_t K = 0; K < wake_pos.size(); K++) {
    k3[K] = Velocity(kernel, pos3[K], wing_pos, wing_gamma, pos3, wake_gamma,
                     forward_flight);
  }
  for (std::size_t n = 0; n < k3.num(); n++) {
    // for trailing edge
    std::fill(k3.iterator_at(n, 0, 0), k3.iterator_at(n, 1, 0),
              -forward_flight);
  }

  // 4th step
  for (std::size_t K = 0; K < wake_pos.size(); K++) {
    pos4[K] = wake_pos[K] + k3[K] * dt;
  }
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (std::size_t K = 0; K < wake_pos.size(); K++) {
    k4[K] = Velocity(kernel, pos4[K], wing_pos, wing_gamma, pos4, wake_gamma,
                     forward_flight);
  }
  for (std::size_t n = 0; n < k4.num(); n++) {
    // for trailing edge
    std::fill(k4.iterator_at(n, 0, 0), k4.iterator_at(n, 1, 0),
              -forward_flight);
  }

  // update
  for (const auto& index : wake_pos.list_index()) {
    std::size_t K, i;
    std::tie(K, std::ignore, i, std::ignore) = index;
    if (i > ignore_steps_) {
      // RK4
      (*next)[K] =
          wake_pos[K] + (k1[K] + k2[K] * 2 + k3[K] * 2 + k3[K]) * dt / 6;
    } else {
      // Euler
      (*next)[K] = wake_pos[K] + k1[K] * dt;
    }
  }
}

void AdamsBashforth2::Advect(MultipleSheet<Eigen::Vector3d>* next,
                             const vortex_kernel::VortexKernel& kernel,
                             const MultipleSheet<Eigen::Vector3d>& wing_pos,
                             const MultipleSheet<double>& wing_gamma,
                             const MultipleSheet<Eigen::Vector3d>& wake_pos,
                             const MultipleSheet<double>& wake_gamma,
                             const Eigen::Vector3d& forward_flight,
                             const double dt) const {
  Velocity(kernel, &v, wing_pos, wing_gamma, wake_pos, wake_gamma, forward_flight);

  if (v_prev.size() == 0) {
    // Euler scheme for the first step
    for (std::size_t K = 0; K < wake_pos.size(); K++) {
      (*next)[K] = wake_pos[K] + v[K] * dt;
    }
  } else {
    // Adams Bashforth 2nd order
    const auto indices = wake_pos.list_index();
    for (const auto& index : indices) {
      std::size_t K, n, i, j;
      std::tie(K, n, i, j) = index;
      if (i >= 1) {
        // Because vortices at trailing edge are shed into wake every step,
        // (i, j) at t corresponds to (i-1, j) at t-dt
        next->at(n, i, j) =
            wake_pos.at(n, i, j) +
            (v.at(n, i, j) * 3 - v_prev.at(n, i - 1, j)) * dt / 2;
      } else {
        // for trailing edge
        next->at(n, i, j) += v[K] * dt;
      }
    }
  }
  // TODO swap
  v_prev.resize(v);
  std::copy(v.begin(), v.end(), v_prev.begin());
}

}  // namespace advect
}  // UVLM

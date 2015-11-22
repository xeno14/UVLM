/**
 * @file simple_simulator.h
 * @brief Add description here
 */
#pragma once

#include "../../proto/uvlm.pb.h"
#include "../morphing.h"
#include "../multiple_sheet/multiple_sheet.h"
#include "../calc_load/joukowski.h"

using multiple_sheet::MultipleSheet;

namespace UVLM {
namespace simulator {

std::vector<Eigen::Vector3d> CollocationPoints(
    const MultipleSheet<Eigen::Vector3d>& pos);

std::vector<Eigen::Vector3d> Normals(const MultipleSheet<Eigen::Vector3d>& pos);

class SimpleSimulator {
 public:
  SimpleSimulator() : forward_flight_(Eigen::Vector3d::Zero()) {}

  /**
   * @brief Add wing information
   *
   * @param chor
   * @todo specify wing generator
   */
  void AddWing(const Morphing& morphing, const double chord, const double span,
               const std::size_t rows, const std::size_t cols,
               const Eigen::Vector3d& origin);

  void Run(const std::size_t steps, const double dt);

  void set_forward_flight(const Eigen::Vector3d& v) { forward_flight_ = v; }

 private:
  MultipleSheet<Eigen::Vector3d> wing_pos;
  MultipleSheet<Eigen::Vector3d> wing_pos_init;
  MultipleSheet<Eigen::Vector3d> wake_pos;
  MultipleSheet<double> wing_gamma;
  MultipleSheet<double> wing_gamma_prev;
  MultipleSheet<double> wake_gamma;
  Eigen::Vector3d forward_flight_;
  std::vector<Morphing> morphings_;

  template <class Range1, class Range2>
  Eigen::MatrixXd CalcMatrix(const std::vector<Eigen::Vector3d>& cpos,
                             const std::vector<Eigen::Vector3d>& normal) const;
  Eigen::VectorXd CalcRhs(const std::vector<Eigen::Vector3d>& cpos,
                          const std::vector<Eigen::Vector3d>& normal,
                          const double t) const;
  Eigen::Vector3d BoundVelocity(const Eigen::Vector3d& x) const;
  Eigen::Vector3d WakeVelocity(const Eigen::Vector3d& x) const;
  Eigen::Vector3d Velocity(const Eigen::Vector3d& x) const {
    return -forward_flight_ + BoundVelocity(x) + WakeVelocity(x);
  }

  void MainLoop(const std::size_t step, const double dt);
};
}  // namespace simulator
}  // namespace UVLM

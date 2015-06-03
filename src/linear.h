
/**
 * @file linear.h
 * @brief Add description here
 */
#pragma once

#include "morphing.h"
#include "uvlm_vortex_ring.h"

#include <vector>

namespace UVLM {
namespace internal {

/**
 * @brief bounded vortex rings から係数行列を求める
 */
Eigen::MatrixXd CalcMatrix(const std::vector<VortexRing>& vortices);

/**
 * @brief 右辺値の主流の寄与分
 */
Eigen::VectorXd CalcRhsUpStream(const Eigen::Vector3d& Vinfty,
                                const std::vector<VortexRing>& vortices);

/**
 * @brief 右辺値のwakeの寄与分
 *
 * 翼上の各渦輪に対し、wakeがつくった流れの法線成分
 */
Eigen::VectorXd CalcRhsWake(const std::vector<VortexRing>& vortices,
                            const std::vector<VortexRing>& wake);

/**
 * @brief 右辺値の変形による寄与分
 */
Eigen::VectorXd CalcRhsMorphing(const std::vector<VortexRing>& vortices,
                                const Morphing& morphing, double t);

}  // namespace internal

inline Eigen::VectorXd SolveLinearProblem(
    const UVLMVortexRing& rings,
    const Eigen::Vector3d Vinfty, const Morphing& morphing, const double t) {
  Eigen::MatrixXd A = internal::CalcMatrix(rings.bound_vortices());
  Eigen::VectorXd rhs =
      internal::CalcRhsUpStream(Vinfty, rings.bound_vortices()) +
      internal::CalcRhsWake(rings.bound_vortices(), rings.wake_vortices()) -
      internal::CalcRhsMorphing(rings.bound_vortices(), morphing, t);
  Eigen::FullPivLU<Eigen::MatrixXd> solver(A);
  return solver.solve(rhs);
}

}  // namespace UVLM

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
template <class InputIterator>
Eigen::MatrixXd CalcMatrix(InputIterator bound_first, InputIterator bound_last);

/**
 * @brief bounded vortex rings から係数行列を求める
 */
inline Eigen::MatrixXd CalcMatrix(const std::vector<VortexRing>& vortices) {
  return CalcMatrix(vortices.cbegin(), vortices.cend());
}

/**
 * @brief 右辺値の主流の寄与分
 */
template <class InputIterator>
Eigen::VectorXd CalcRhsUpStream(const Eigen::Vector3d& Vinfty,
                                InputIterator bound_first,
                                InputIterator bound_last);

/**
 * @brief 右辺値の主流の寄与分
 */
inline Eigen::VectorXd CalcRhsUpStream(const Eigen::Vector3d& Vinfty,
                                const std::vector<VortexRing>& vortices) {
  return CalcRhsUpStream(Vinfty, vortices.cbegin(), vortices.cend());
}

/**
 * @brief 右辺値のwakeの寄与分
 *
 * 翼上の各渦輪に対し、wakeがつくった流れの法線成分
 */
template <class InputIterator1, class InputIterator2>
Eigen::VectorXd CalcRhsWake(InputIterator1 bound_first,
                            InputIterator1 bound_last,
                            InputIterator2 wake_first,
                            InputIterator2 wake_last);

/**
 * @brief 右辺値のwakeの寄与分
 *
 * 翼上の各渦輪に対し、wakeがつくった流れの法線成分
 */
inline Eigen::VectorXd CalcRhsWake(const std::vector<VortexRing>& vortices,
                            const std::vector<VortexRing>& wake) {
  return CalcRhsWake(vortices.cbegin(), vortices.cend(), wake.cbegin(),
                     wake.cend());
}

/**
 * @brief 右辺値の変形による寄与分
 */
template <class InputIterator1>
Eigen::VectorXd CalcRhsMorphing(InputIterator1 bound_first,
                                InputIterator1 bound_last,
                                const Morphing& morphing, double t);

/**
 * @brief 右辺値の変形による寄与分
 */
inline Eigen::VectorXd CalcRhsMorphing(const std::vector<VortexRing>& vortices,
                                const Morphing& morphing, double t) {
  return CalcRhsMorphing(vortices.cbegin(), vortices.cend(), morphing, t);
}

}  // namespace internal

/** @brief 連立方程式を解いて循環を求める
 */
template <class InputIterator1, class InputIterator2>
inline Eigen::VectorXd SolveLinearProblem(
    InputIterator1 bound_first, InputIterator1 bound_last,
    InputIterator2 wake_first, InputIterator2 wake_last,
    const Eigen::Vector3d Vinfty, const Morphing& morphing, const double t) {
  Eigen::MatrixXd A = internal::CalcMatrix(bound_first, bound_last);
  Eigen::VectorXd rhs = 
      internal::CalcRhsUpStream(Vinfty, bound_first, bound_last) +
      internal::CalcRhsWake(bound_first, bound_last, wake_first, wake_last) -
      internal::CalcRhsMorphing(bound_first, bound_last, morphing, t);
  Eigen::FullPivLU<Eigen::MatrixXd> solver(A);
  return solver.solve(rhs);
}

}  // namespace UVLM

#include "linear.inl"

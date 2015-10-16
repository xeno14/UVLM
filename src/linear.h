/**
 * @file linear.h
 * @brief Add description here
 */
#pragma once

#include "morphing.h"
#include "vortex_container.h"

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
 * @brief Calc velocity of morphing at each collocation point on panels
 */
template <class InputIterator, class OutputIterator>
OutputIterator CalcMorphingVelocityOnWing(InputIterator panel_first,
                                          InputIterator panel_last,
                                          const Morphing& morphing, double t,
                                          OutputIterator result);

/**
 * @brief 右辺値の変形による寄与分
 */
inline Eigen::VectorXd CalcRhsMorphing(const std::vector<VortexRing>& vortices,
                                const Morphing& morphing, double t) {
  return CalcRhsMorphing(vortices.cbegin(), vortices.cend(), morphing, t);
}

/**
 * @brief Calc right hand contributed from morphings 
 * @pre containers and morphings have exactly same size
 */
inline Eigen::VectorXd CalcRhsMorphings(
    const std::vector<VortexContainer>& containers,
    const std::vector<Morphing>& morphings, double t);

}  // namespace internal

/** @brief 連立方程式を解いて循環を求める
 */
template <class InputIterator1, class InputIterator2>
Eigen::VectorXd SolveLinearProblem(InputIterator1 bound_first,
                                   InputIterator1 bound_last,
                                   InputIterator2 wake_first,
                                   InputIterator2 wake_last,
                                   const Eigen::Vector3d Vinfty,
                                   const Morphing& morphing, const double t);

/**
 * @brief Solve simultaneous equations obtained from no-penetration condition.
 */
Eigen::VectorXd SolveLinearProblem(
    const std::vector<VortexContainer>& containers,
    const std::vector<Morphing>& morphings, const Eigen::Vector3d& freestream,
    const double t);

}  // namespace UVLM

#include "linear.inl"

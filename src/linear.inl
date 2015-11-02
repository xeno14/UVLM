/**
 * @file linear.cpp
 * @brief Add description here
 */
#pragma once

#include "shed.h"

#include <glog/logging.h>

namespace UVLM {
namespace internal {

template <class InputIterator>
Eigen::MatrixXd CalcMatrix(InputIterator bound_first,
                           InputIterator bound_last) {
  const std::size_t len = std::distance(bound_first, bound_last);
  Eigen::MatrixXd res(len, len);

  Eigen::Vector3d velocity;
  for (std::size_t i = 0; i < len; i++) {
    for (std::size_t j = 0; j < len; j++) {
      const auto& vi = *(bound_first + i);
      const auto& vj = *(bound_first + j);
      vj.BiotSavartLaw(&velocity,
                       vi.Centroid(),  // 渦iの位置
                       1.);
      res(i, j) = velocity.dot(vi.Normal());  // 渦jが渦iに及ぼす影響
    }
  }
  return res;
}

template <class InputIterator>
Eigen::VectorXd CalcRhsFreeStream(const Eigen::Vector3d& freestream,
    InputIterator bound_first,
    InputIterator bound_last) {
  const std::size_t len = std::distance(bound_first, bound_last);
  Eigen::VectorXd res(len);
  for (std::size_t i = 0; i < len; i++) {
    const auto& v = *(bound_first + i);
    res(i) = freestream.dot(v.Normal());
  }
  return res;
}

template <class InputIterator1, class InputIterator2>
Eigen::VectorXd CalcRhsWake(InputIterator1 bound_first,
                            InputIterator1 bound_last,
                            InputIterator2 wake_first,
                            InputIterator2 wake_last) {
  const std::size_t len = std::distance(bound_first, bound_last);
  Eigen::VectorXd res(len);
  Eigen::Vector3d vel;
  for (std::size_t i = 0; i < len; i++) {
    const auto& v = *(bound_first + i);
    InducedVelocity(&vel, v.Centroid(), wake_first, wake_last);
    res(i) = vel.dot(v.Normal());
  }
  return res;
}

template <class InputIterator1>
Eigen::VectorXd CalcRhsMorphing(InputIterator1 bound_first,
                                InputIterator1 bound_last,
                                const Morphing& morphing, double t) {
  const std::size_t len = std::distance(bound_first, bound_last);
  Eigen::VectorXd res(len);
  Eigen::Vector3d vel;
  for (std::size_t i = 0; i < len; i++) {
    const auto& v = *(bound_first + i);
    morphing.Velocity(&vel, v.ReferenceCentroid(), t);
    res(i) = vel.dot(v.Normal());
  }
  return res;
}

template <class InputIterator, class OutputIterator>
OutputIterator CalcMorphingVelocityOnWing(InputIterator panel_first,
                                          InputIterator panel_last,
                                          const Morphing& morphing, double t,
                                          OutputIterator result) {
  Eigen::Vector3d Vls;  // velocity of lifting surface

  for (auto panel = panel_first; panel != panel_last; ++panel, ++result) {
    morphing.Velocity(&Vls, panel->ReferenceCentroid(), t);
    *result = Vls.dot(panel->Normal());
  }
  return result;
}

inline Eigen::VectorXd CalcRhsMorphings(
    const std::vector<VortexContainer>& containers,
    const std::vector<Morphing>& morphings, double t) {
  const std::size_t wake_offset = 
      ::UVLM::CountTotalSize(containers.begin(), containers.end());
  if (containers.size() != morphings.size()) {
    LOG(FATAL) << "containers and morphings have to be same size";
  }
  Eigen::VectorXd res(wake_offset);
  std::vector<double> rhs(wake_offset);
  auto it = rhs.begin();
  for (std::size_t i=0; i<containers.size(); i++) {
    it = CalcMorphingVelocityOnWing(containers[i].cbegin(),
                                    containers[i].cend(), morphings[i], t, it);
  }
  std::copy(rhs.begin(), rhs.end(), res.data());
  return res;
}

}  // namespace internal

template <class InputIterator1, class InputIterator2>
Eigen::VectorXd SolveLinearProblem(InputIterator1 bound_first,
                                   InputIterator1 bound_last,
                                   InputIterator2 wake_first,
                                   InputIterator2 wake_last,
                                   const Eigen::Vector3d freestream,
                                   const Morphing& morphing, const double t) {
  Eigen::MatrixXd A = internal::CalcMatrix(bound_first, bound_last);
  Eigen::VectorXd rhs =
      internal::CalcRhsMorphing(bound_first, bound_last, morphing, t)
      - internal::CalcRhsFreeStream(freestream, bound_first, bound_last);
      - internal::CalcRhsWake(bound_first, bound_last, wake_first, wake_last);
  Eigen::FullPivLU<Eigen::MatrixXd> solver(A);
  if (!solver.isInvertible()) {
    LOG(FATAL) << "Matrix is not invertible.";
  }
  return solver.solve(rhs);
}

Eigen::VectorXd SolveLinearProblem(
    const std::vector<VortexContainer>& containers,
    const std::vector<Morphing>& morphings, const Eigen::Vector3d& freestream,
    const double t) {
  // All vortex rings
  const auto& vortices = containers.begin()->vortices();
  const std::size_t wake_offset =
      ::UVLM::CountTotalSize(containers.begin(), containers.end());

  auto bound_first = vortices->begin();
  auto bound_last = vortices->begin() + wake_offset;
  auto wake_first = vortices->begin() + wake_offset;
  auto wake_last = vortices->end();

  Eigen::MatrixXd A =
      internal::CalcMatrix(bound_first, bound_last);
  Eigen::VectorXd rhs =
      internal::CalcRhsMorphings(containers, morphings, t)
      - internal::CalcRhsFreeStream(freestream, bound_first, bound_last);
      - internal::CalcRhsWake(bound_first, bound_last, wake_first, wake_last);
  Eigen::FullPivLU<Eigen::MatrixXd> solver(A);
  if (!solver.isInvertible()) {
    LOG(FATAL) << "Matrix is not invertible.";
  }
  return solver.solve(rhs);
}

}  // namespace UVLM

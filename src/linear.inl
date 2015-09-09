/**
 * @file linear.cpp
 * @brief Add description here
 */
#pragma once


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
Eigen::VectorXd CalcRhsUpStream(const Eigen::Vector3d& Vinfty,
    InputIterator bound_first,
    InputIterator bound_last) {
  const std::size_t len = std::distance(bound_first, bound_last);
  Eigen::VectorXd res(len);
  for (std::size_t i = 0; i < len; i++) {
    res(i) = Vinfty.dot((bound_first + i)->Normal());
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
    res(i) = 0;
    const auto& v = *(bound_first + i);
    for (auto wake = wake_first; wake != wake_last; ++wake) {
      wake->BiotSavartLaw(&vel, v.Centroid());  // 渦輪の中心の流速
      res(i) += vel.dot(v.Normal());
    }
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

}  // namespace internal

template <class InputIterator1, class InputIterator2>
Eigen::VectorXd SolveLinearProblem(InputIterator1 bound_first,
                                   InputIterator1 bound_last,
                                   InputIterator2 wake_first,
                                   InputIterator2 wake_last,
                                   const Eigen::Vector3d Vinfty,
                                   const Morphing& morphing, const double t) {
  Eigen::MatrixXd A = internal::CalcMatrix(bound_first, bound_last);
  Eigen::VectorXd rhs =
      internal::CalcRhsUpStream(Vinfty, bound_first, bound_last) +
      internal::CalcRhsWake(bound_first, bound_last, wake_first, wake_last) -
      internal::CalcRhsMorphing(bound_first, bound_last, morphing, t);
  Eigen::PartialPivLU<Eigen::MatrixXd> solver(A);
  return solver.solve(rhs);
}

}  // namespace UVLM

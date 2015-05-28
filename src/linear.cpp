/**
 * @file linear.cpp
 * @brief Add description here
 */

#include "linear.h"


namespace UVLM {
namespace internal {

Eigen::MatrixXd CalcMatrix(const std::vector<VortexRing>& vortices) {
  Eigen::MatrixXd res(vortices.size(), vortices.size());

  Eigen::Vector3d tmp;
  for (std::size_t i=0; i<vortices.size(); i++) {
    for (std::size_t j=0; j<vortices.size(); j++) {
      const auto& vi = vortices[i];
      const auto& vj = vortices[j];
      vj.BiotSavartLaw(
        &tmp,
        vi.Centroid(),    // 渦iの位置
        1.);
      res(i, j) = tmp.dot(vi.Normal());  // 渦jが渦iに及ぼす影響
    }
  }
  return res;
}

Eigen::VectorXd CalcRhsUpStream(const Eigen::Vector3d& Vinfty,
                                const std::vector<VortexRing>& vortices) {
  Eigen::VectorXd res(vortices.size());
  for (std::size_t i = 0; i < vortices.size(); i++) {
    res(i) = Vinfty.dot(vortices[i].Normal());
  }
  return res;
}

Eigen::VectorXd CalcRhsWake(const std::vector<VortexRing>& vortices,
                            const std::vector<VortexRing>& wake) {
  Eigen::VectorXd res(vortices.size());
  Eigen::Vector3d vel;
  for (std::size_t i = 0; i < vortices.size(); i++) {
    res(i) = 0;
    const auto& v = vortices[i];
    for (const auto& w : wake) {
      w.BiotSavartLaw(&vel, v.Centroid());  // 渦輪の中心の流速
      res(i) += vel.dot(v.Normal());
    }
  }
  return res;
}

}  // namespace internal

}  // namespace UVLM

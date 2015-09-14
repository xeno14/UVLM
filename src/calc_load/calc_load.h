
/**
 * @file calc_load.h
 * @brief Add description here
 */
#pragma once

#include "../vortex_container.h"
#include "../shed.h"

namespace UVLM {

/**
 * パネルi,jでの圧力差を求める
 */
template <class InputIterator>
double CalcDP(const std::size_t i, const std::size_t j,
              const VortexContainer& vb, const VortexContainer& vb_prev,
              InputIterator wake_first, InputIterator wake_last,
              const Eigen::Vector3d& Vinfty, const double rho,
              const double dt) {
  Eigen::Vector3d vw;  // wakeによって作られた流れ
  InducedVelocity(&vw, vb.at(i, j).Centroid(), wake_first, wake_last);
  const Eigen::Vector3d V = Vinfty + vw;

  // 端っこのパネルでは差を使わない
  const double dg_dx =
      (i == 0 ? vb.at(i, j).gamma()
              : vb.at(i, j).gamma() - vb.at(i - 1, j).gamma()) /
      vb.at(i, j).CalcC();
  const double dg_dy =
      (j == 0 ? vb.at(i, j).gamma()
              : vb.at(i, j).gamma() - vb.at(i, j - 1).gamma()) /
      vb.at(i, j).CalcB();
  const double dg_dt = (vb.at(i, j).gamma() - vb_prev.at(i, j).gamma()) / dt;

  return rho * (V.dot(vb.at(i, j).TanVecChord()) * dg_dx +
                V.dot(vb.at(i, j).TanVecSpan()) * dg_dy + dg_dt);
}

template <class InputIterator>
Eigen::Vector3d CalcLoad(const VortexContainer& vb, const VortexContainer& vb_prev,
                InputIterator wake_first, InputIterator wake_last,
                const Eigen::Vector3d& Vinfty, const double rho,
                const double dt) {
  Eigen::Vector3d res = Eigen::Vector3d::Zero();
  for (std::size_t i=0; i<vb.rows(); i++) {
    for (std::size_t j=0; j<vb.cols(); j++) {
      const double deltaP = CalcDP(i, j, vb, vb_prev, wake_first, wake_last, Vinfty, rho, dt);
      const Eigen::Vector3d n = vb.at(i, j).Normal();
      const double dS = vb.at(i, j).CalcB() * vb.at(i, j).CalcC();

      res += n * deltaP* dS;
    }
  }
  return res;
}

}  // namespace UVLM

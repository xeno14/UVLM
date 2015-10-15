
/**
 * @file joukowski.h
 * @brief Add description here
 */
#pragma once

#include "../../proto/uvlm.pb.h"
#include "../vortex_container.h"
#include "../shed.h"
#include "../util.h"
#include "../morphing.h"
#include "calc_load.h"

namespace UVLM {
namespace calc_load {
namespace internal {

// unstready part
template <class InputIterator>
inline void JoukowskiSteadyOnPanel(Eigen::Vector3d* result, const VortexRing& v,
                                   InputIterator vortices_first,
                                   InputIterator vortices_last,
                                   const Morphing& morphing,
                                   const Eigen::Vector3d& freestream,
                                   const double rho, const double t) {
  Eigen::Vector3d mid, mid0;
  Eigen::Vector3d U;     // velocity on the point
  Eigen::Vector3d Um;    // velocity contribution from the surface motion
  Eigen::Vector3d Uind;  // induced velocity from all vortices
  *result = Eigen::Vector3d::Zero();
  v.ForEachSegment([&](const auto& start, const auto& end, const auto& start0,
                       const auto& end0) {
    mid = (start + end) / 2;
    mid0 = (start0 + end0) / 2;
    UVLM::InducedVelocity(&Uind, mid, vortices_first, vortices_last);
    // U = (Ub + Uw) + Um = Uind + Um
    Um = internal::CalcUm(morphing, mid0, freestream, t);
    U = Uind + Um;
    *result += U.cross(end - start) * rho * v.gamma();
  });
}

// unstready part
inline void JoukowskiUnsteadyOnPanel(Eigen::Vector3d* result,
                                     const VortexRing& v,
                                     const VortexRing& v_prev, const double rho,
                                     const double dt) {
  const double dGdt = (v.gamma() - v_prev.gamma()) / dt;
  *result = v.Normal() * rho * dGdt * v.CalcB() * v.CalcC();
}

}  // namespace internal

template <class InputIterator>
AerodynamicLoad CalcLoadJoukowski(
    const VortexContainer& vb, const VortexContainer& vb_prev,
    InputIterator wake_first, InputIterator wake_last, const Morphing& morphing,
    const Eigen::Vector3d& freestream, const double rho, const double t,
    const double dt) {
  const auto& vortices = *vb.vortices();
  auto dim = DoubleLoop(vb.rows(), vb.cols());
  double Fx=0, Fy=0, Fz=0;
  std::size_t index;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:Fx, Fy, Fz)
#endif
  for (index = 0; index < dim.size(); ++index) {
    const auto i = dim[index].first;
    const auto j = dim[index].second;
    Eigen::Vector3d dF_st = Eigen::Vector3d::Zero();
    Eigen::Vector3d dF_unst = Eigen::Vector3d::Zero();
    internal::JoukowskiSteadyOnPanel(&dF_st, vb.at(i, j), vortices.cbegin(),
                                     vortices.cend(), morphing, freestream, rho,
                                     t);
    internal::JoukowskiUnsteadyOnPanel(&dF_unst, vb.at(i, j), vb_prev.at(i, j),
                                       rho, dt);
    Fx += dF_st.x() + dF_unst.x();
    Fy += dF_st.y() + dF_unst.y();
    Fz += dF_st.z() + dF_unst.z();
  }
  Eigen::Vector3d F(Fx, Fy, Fz);
  return AerodynamicLoad{F, 0, 0};
}

}  // namespace calc_load
}  // namespace UVLM

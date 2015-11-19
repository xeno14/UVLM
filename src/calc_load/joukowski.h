
/**
 * @file joukowski.h
 * @brief Add description here
 */
#pragma once

#include "../../proto/uvlm.pb.h"
#include "../multiple_sheet/multiple_sheet.h"
#include "../vortex_container.h"
#include "../shed.h"
#include "../util.h"
#include "../morphing.h"
#include "functions.h"

using multiple_sheet::MultipleSheet;

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
                                   const double rho, const double t,
                                   bool edge_flag) {
  Eigen::Vector3d mid, mid0;
  Eigen::Vector3d U;     // velocity on the point
  Eigen::Vector3d Um;    // velocity contribution from the surface motion
  Eigen::Vector3d Uind;  // induced velocity from all vortices
  *result = Eigen::Vector3d::Zero();
  for (std::size_t i = 0; i < v.nodes().size(); i++) {
    if (edge_flag && i == 1) continue;  // ignore trailing edge_end
    const Vector3d& start = v.nodes()[(i + 1) % v.nodes().size()];
    const Vector3d& end = v.nodes()[i];
    const Vector3d& start0 = v.nodes0()[(i + 1) % v.nodes().size()];
    const Vector3d& end0 = v.nodes0()[i];
    mid = (start + end) / 2;
    mid0 = (start0 + end0) / 2;
    UVLM::InducedVelocity(&Uind, mid, vortices_first, vortices_last);
    // U = (Ub + Uw) + Um = Uind + Um
    Um = internal::CalcUm(morphing, mid0, freestream, t);
    U = Uind + Um;
    *result += U.cross(end - start) * rho * v.gamma();
  };
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

inline AerodynamicLoad CalcLoadJoukowski(const VortexContainer& vb,
                                  const VortexContainer& vb_prev,
                                  const Morphing& morphing,
                                  const Eigen::Vector3d& freestream,
                                  const double rho, const double t,
                                  const double dt) {
  const auto& vortices = *vb.vortices();
  auto dim = DoubleLoop(vb.rows(), vb.cols());
  double Fx=0, Fy=0, Fz=0, Pin=0, Pout=0;

#ifdef _OPENMP
#pragma omp parallel for reduction(+:Fx, Fy, Fz, Pin)
#endif
  for (std::size_t index = 0; index < dim.size(); ++index) {
    const auto i = dim[index].first;
    const auto j = dim[index].second;
    Eigen::Vector3d dF_st = Eigen::Vector3d::Zero();
    Eigen::Vector3d dF_unst = Eigen::Vector3d::Zero();


    internal::JoukowskiSteadyOnPanel(&dF_st, vb.at(i, j), vortices.cbegin(),
                                     vortices.cend(), morphing, freestream, rho,
                                     t, i+1==vb.rows());
    internal::JoukowskiUnsteadyOnPanel(&dF_unst, vb.at(i, j), vb_prev.at(i, j),
                                       rho, dt);
    const Eigen::Vector3d dF = dF_st + dF_unst;
    Fx += dF.x();
    Fy += dF.y();
    Fz += dF.z();

      // auto pos = vb.at(i,j).Centroid();
      // std::cout << pos.x() << " " << pos.y() << " " << pos.z() << " " << dF.x()
      //   << " " << dF.y() << " " << dF.z() << " " << vb.at(i,j).gamma() << std::endl;

    // TODO duplicate Joukowski
    Eigen::Vector3d Vls;    // velocity of motion
    morphing.Velocity(&Vls, vb.at(i, j).ReferenceCentroid(), t);
    Pin += dF.dot(Vls);
  }
  // std::cout << "\n\n";
  Eigen::Vector3d F(Fx, Fy, Fz);
  Pout = Fx * freestream.norm() * (-1);    // assume foward flight
  return AerodynamicLoad{F, Pin, Pout};
}

struct VortexLine {
  Eigen::Vector3d p0, p1;
  Eigen::Vector3d p0_init, p1_init;
  double gamma;
};

std::vector<VortexLine> GetLines(const MultipleSheet<Eigen::Vector3d>& pos,
                                 const MultipleSheet<Eigen::Vector3d>& pos_init,
                                 const MultipleSheet<double>& gamma,
                                 std::size_t n);

}  // namespace calc_load
}  // namespace UVLM

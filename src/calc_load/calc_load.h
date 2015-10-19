
/**
 * @file calc_load.h
 * @brief Calculation of aerodynamic loads.
 *
 * Formulations are followed "Induced-Drag Calculation in the Unsteady Vortex
 * Lattice Method" (2013).
 */
#pragma once

#include "../../proto/uvlm.pb.h"
#include "../vortex_container.h"
#include "../shed.h"
#include "../morphing.h"

namespace UVLM {
namespace calc_load {
namespace internal {

inline Eigen::Matrix3d CalcProjectionOperator(const Eigen::Vector3d& Um) {
  Eigen::Vector3d Um_ = Um; // normalized vector
  Um_.normalize();
  return Eigen::Matrix3d::Identity() - Um_ * Um_.transpose();
}

/**
 * @brief Calc velocity contribution from the surface motion
 */
inline Eigen::Vector3d CalcUm(const Morphing& morphing,
                              const Eigen::Vector3d& x0, 
                              const Eigen::Vector3d& freestream,
                              const double t) {
  Eigen::Vector3d Uls;
  morphing.Velocity(&Uls, x0, t);
  return freestream - Uls;
}

inline double CalcChordwiseDGamma(const std::size_t i, const std::size_t j,
                                  const VortexContainer& vb) {
  return i == 0 ? vb.at(i, j).gamma()
                : vb.at(i, j).gamma() - vb.at(i - 1, j).gamma();
}

inline double CalcSpanwiseDGamma(const std::size_t i, const std::size_t j,
                                 const VortexContainer& vb) {
  return j == 0 ? vb.at(i, j).gamma()
                : vb.at(i, j).gamma() - vb.at(i, j - 1).gamma();
}

inline double CalcGammaTimeDerivate(const VortexRing& v, const VortexRing& v_prev,
    const double dt) {
  return (v.gamma() - v_prev.gamma()) / dt;
}

inline double CalcLocalLift(const std::size_t i, const std::size_t j,
                           const VortexContainer& vb,
                           const VortexContainer& vb_prev,
                           const Eigen::Vector3d& Um, const Eigen::Vector3d& Uw,
                           const double alpha,
                           const double rho, const double dt) {
  const VortexRing& v = vb.at(i, j);
  const VortexRing& v_prev = vb_prev.at(i, j);
  const double b = v.CalcB();
  const double c = v.CalcC();
  // -1: direction of gamma is opposite to the articles
  return rho * b * c * cos(alpha) *
         ((Um + Uw).dot(v.TanVecChord()) * CalcChordwiseDGamma(i, j, vb) / c +
          (Um + Uw).dot(v.TanVecSpan()) * CalcSpanwiseDGamma(i, j, vb) / b +
          CalcGammaTimeDerivate(v, v_prev, dt));
}

inline double CalcLocalDrag(const std::size_t i, const std::size_t j,
                           const VortexContainer& vb,
                           const VortexContainer& vb_prev,
                           const Eigen::Matrix3d& P,
                           const Eigen::Vector3d& Ubc, const Eigen::Vector3d& Uw,
                           const double alpha,
                           const double rho, const double dt) {
  const VortexRing& v = vb.at(i, j);
  const VortexRing& v_prev = vb_prev.at(i, j);
  const double b = v.CalcB();
  const double c = v.CalcC();
  return rho * (
      -(Ubc + Uw).dot(P * v.Normal()) * CalcChordwiseDGamma(i, j, vb) * b +
      CalcGammaTimeDerivate(v, v_prev, dt) * b * c * sin(alpha));
}

}  // namespace internal

struct AerodynamicLoad {
  Eigen::Vector3d F;
  double Pin, Pout;
};

template <class InputIterator>
AerodynamicLoad CalcLoad(const VortexContainer& vb,
                         const VortexContainer& vb_prev,
                         InputIterator wake_first, InputIterator wake_last,
                         const Morphing& morphing, const Eigen::Vector3d& freestream,
                         const double rho, const double t, const double dt) {
  Eigen::Vector3d F = Eigen::Vector3d::Zero();
  double Pin = 0;
  for (std::size_t i = 0; i < vb.rows(); i++) {
    for (std::size_t j = 0; j < vb.cols(); j++) {
      const Eigen::Vector3d centroid = vb.at(i, j).Centroid();
      const Eigen::Vector3d normal = vb.at(i, j).Normal();

      const Eigen::Vector3d Um =
          internal::CalcUm(morphing, centroid, freestream, t);
      const Eigen::Matrix3d P = internal::CalcProjectionOperator(Um);
      Eigen::Vector3d Ubc, Uw;
      ChordwiseInducedVelocity(&Ubc, centroid, vb.cbegin(), vb.cend());
      InducedVelocity(&Uw, centroid, wake_first, wake_last);
      const double alpha = vb.at(i, j).AngleOfAttack(Um);

      const double Llocal =
          internal::CalcLocalLift(i, j, vb, vb_prev, Um, Uw, alpha, rho, dt);
      const double Dlocal = internal::CalcLocalDrag(i, j, vb, vb_prev, P, Ubc,
                                                    Uw, alpha, rho, dt);

      Eigen::Vector3d Um_ = Um;
      Um_.normalize();

      Eigen::Vector3d dF = Um_ * Dlocal + P * normal * Llocal;
      F += dF;
      Pin += dF.dot(Um) * (-1);
    }
  }
  return AerodynamicLoad {F, Pin, F.x() * freestream.norm() * (-1)};
}

}  // namespace calc_load

namespace calc_load {
namespace internal {

proto::Snapshot2 ReadSnapshot2(const std::string& filename);
std::vector<proto::Snapshot2> ReadSnapshot2(
    const std::vector<std::string>& filenames);
std::vector<Eigen::Vector3d> Calc(const proto::Snapshot2& s0, const proto::Snapshot2& s1);

}  // namespace internal

std::vector<std::string> Snapshot2Paths(const std::string& pattern);

void Start(const ::UVLM::Morphing& m, const std::string& output_path);

}  // namespace calc_load
}  // namespace UVLM

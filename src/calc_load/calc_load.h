
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
  Eigen::Vector3d res;
  morphing.Velocity(&res, x0, t);
  res = freestream - res;
  return res;
}

}  // namespace internal
}  // namespace calc_load

inline void MatrixForLocalUnitVector(Eigen::Matrix3d* m,
    const Eigen::Vector3d& n, const Eigen::Vector3d& t, const double alpha) {
  Eigen::Vector3d q = t.cross(n);
  Eigen::Matrix3d T;
  T << t.x(), n.x(), q.x(),
       t.y(), n.y(), q.y(),
       t.z(), n.z(), q.z();
  // In Ghommem's article, T_alpha here is transposed.
  Eigen::Matrix3d T_alpha;
  T_alpha << cos(alpha), -sin(alpha), 0,
             sin(alpha),  cos(alpha), 0,
              0,          0,          1;
  *m = T * T_alpha * T.transpose();
}

inline void LocalUnitVector(Eigen::Vector3d* e_lift, Eigen::Vector3d* e_drag,
                            const Eigen::Vector3d& n, const Eigen::Vector3d& t,
                            const double alpha) {
  Eigen::Matrix3d m;
  MatrixForLocalUnitVector(&m, n, t, alpha);
  *e_lift = m * n;
  *e_drag = m * t;
}

template <class InputIterator>
double CalcDragOnPanel(const std::size_t i, const std::size_t j,
              const VortexContainer& vb, const VortexContainer& vb_prev,
              InputIterator wake_first, InputIterator wake_last,
              const Morphing& morphing, const Eigen::Vector3d& inlet,
              const double rho, const double t, const double dt) {
  // calc for angle of attack
  Eigen::Vector3d V_kinematic = calc_load::internal::CalcUm(
      morphing, vb.at(i, j).ReferenceCentroid(), inlet, t);
  const double alpha = vb.at(i, j).AngleOfAttack(V_kinematic);

  // calc for induced velocity due to wake vortices
  Eigen::Vector3d V_wake;
  const Eigen::Vector3d centroid = vb.at(i, j).Centroid();
  InducedVelocity(&V_wake, centroid, wake_first, wake_last);

  // Calc for V_ind
  Eigen::Vector3d V_ind;
  ChordwiseInducedVelocity(&V_ind, centroid, vb.cbegin(), vb.cend());

  const double induced = vb.at(i, j).Normal().dot((V_ind + V_wake));

  const double C = vb.at(i, j).CalcC();
  const double B = vb.at(i, j).CalcB();
  const double dS = C * B;
  const double dg_dx =
      (i == 0 ? vb.at(i, j).gamma()
              : vb.at(i, j).gamma() - vb.at(i - 1, j).gamma()) / C;
  const double dg_dt = (vb.at(i, j).gamma() - vb_prev.at(i, j).gamma()) / dt;

  return rho * dS * (induced * dg_dx + dg_dt * sin(alpha));
}

inline double CalcPinOnPanel(const Eigen::Vector3d& dF,
                             const Eigen::Vector3d& normal,
                             const Eigen::Vector3d& V_motion) {
  return dF.dot(normal) * normal.dot(V_motion);
}

inline double CalcPout(const Eigen::Vector3d& F, const Eigen::Vector3d& inlet) {
  Eigen::Vector3d unit = inlet;
  unit.normalize();
  return inlet.norm() * F.dot(unit);
}

struct AerodynamicLoad {
  Eigen::Vector3d F;
  double Pin, Pout;
};

template <class InputIterator>
void CalcLoadOnPanel(Eigen::Vector3d* dL, Eigen::Vector3d* dD, double* dPin,
                     const std::size_t i, const std::size_t j,
                     const VortexContainer& vb, const VortexContainer& vb_prev,
                     InputIterator wake_first, InputIterator wake_last,
                     const Morphing& morphing, const Eigen::Vector3d& inlet,
                     const double rho, const double t, const double dt) {
  // calc for angle of attack
  Eigen::Vector3d V_kinematic;
  morphing.Velocity(&V_kinematic, vb.at(i, j).ReferenceCentroid(), t);
  V_kinematic -= inlet;   // -inlet is equivalent to forward flight
  const double alpha = vb.at(i, j).AngleOfAttack(V_kinematic);

  // calc for unit vectors
  const Eigen::Vector3d n = vb.at(i, j).Normal();
  Eigen::Vector3d e_lift, e_drag;
  LocalUnitVector(&e_lift, &e_drag, n, vb.at(i, j).Tangent(), alpha);

  // calc for induced velocity due to wake vortices
  Eigen::Vector3d V_wake;
  const Eigen::Vector3d centroid = vb.at(i, j).Centroid();
  InducedVelocity(&V_wake, centroid, wake_first, wake_last);

  // Calc for V_ind
  Eigen::Vector3d V_ind;
  ChordwiseInducedVelocity(&V_ind, centroid, vb.cbegin(), vb.cend());

  const double induced = n.dot((V_ind + V_wake));
  const double v_k = n.dot(V_kinematic);

  const double C = vb.at(i, j).CalcC();
  const double B = vb.at(i, j).CalcB();
  const double dS = C * B;
  const double dg_dx =
      (i == 0 ? vb.at(i, j).gamma()
              : vb.at(i, j).gamma() - vb.at(i - 1, j).gamma()) / C;
  const double dg_dt = (vb.at(i, j).gamma() - vb_prev.at(i, j).gamma()) / dt;

  *dL = e_lift * rho * dS * (v_k * dg_dx + dg_dt) * cos(alpha);
  *dD = e_drag * rho * dS * (induced * dg_dx + dg_dt * sin(alpha));
  *dPin = CalcPinOnPanel(*dL + *dD, n, V_kinematic);
}

template <class InputIterator>
auto CalcLoad(const VortexContainer& vb,
                         const VortexContainer& vb_prev,
                         InputIterator wake_first, InputIterator wake_last,
                         const Morphing& morphing, const Eigen::Vector3d& inlet,
                         const double rho, const double t, const double dt) {
  Eigen::Vector3d F = Eigen::Vector3d::Zero();
  double Pin=0, Pout;
  for (std::size_t i = 0; i < vb.rows(); i++) {
    for (std::size_t j = 0; j < vb.cols(); j++) {
      Eigen::Vector3d dL, dD;
      double dPin;
      CalcLoadOnPanel(&dL, &dD, &dPin, i, j, vb, vb_prev, wake_first, wake_last,
                      morphing, inlet, rho, t, dt);
      const Eigen::Vector3d dF = dL + dD;
      F += dF;
      Pin += dPin;
    }
  }
  Pout = CalcPout(F, inlet);

  AerodynamicLoad res {F, Pin, Pout};
  return res;
}

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

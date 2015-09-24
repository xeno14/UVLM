
/**
 * @file calc_load.h
 * @brief Add description here
 */
#pragma once

#include "../../proto/uvlm.pb.h"
#include "../vortex_container.h"
#include "../shed.h"
#include "../morphing.h"

namespace UVLM {

inline void MatrixForLocalUnitVector(Eigen::Matrix3d* m,
    const Eigen::Vector3d& n, const Eigen::Vector3d& t, const double alpha) {
  Eigen::Vector3d q = t.cross(n);
  Eigen::Matrix3d T;
  T << t.x(), n.x(), q.x(),
       t.y(), n.y(), q.y(),
       t.z(), n.z(), q.z();
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
  Eigen::Vector3d V_kinematic;
  morphing.Velocity(&V_kinematic, vb.at(i, j).ReferenceCentroid(), t);
  V_kinematic -= inlet;   // -inlet = forward flight
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


/**
 * パネルi,jでの圧力差を求める
 */
template <class InputIterator>
double CalcDP(const std::size_t i, const std::size_t j,
              const VortexContainer& vb, const VortexContainer& vb_prev,
              InputIterator wake_first, InputIterator wake_last,
              const Morphing& morphing, const Eigen::Vector3d& inlet,
              const double rho, const double t, const double dt) {
  Eigen::Vector3d vw;  // wakeによって作られた流れ
  Eigen::Vector3d centroid = vb.at(i, j).Centroid();
  InducedVelocity(&vw, centroid, wake_first, wake_last);
  Eigen::Vector3d v_morphing;
  morphing.Velocity(&v_morphing, centroid, t);
  const Eigen::Vector3d V = inlet - v_morphing + vw;

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
Eigen::Vector3d CalcLoad(const VortexContainer& vb,
                         const VortexContainer& vb_prev,
                         InputIterator wake_first, InputIterator wake_last,
                         const Morphing& morphing, const Eigen::Vector3d& inlet,
                         const double rho, const double t, const double dt) {
  Eigen::Vector3d res = Eigen::Vector3d::Zero();
  double drag = 0;
  for (std::size_t i = 0; i < vb.rows(); i++) {
    for (std::size_t j = 0; j < vb.cols(); j++) {
      const double deltaP = CalcDP(i, j, vb, vb_prev, wake_first, wake_last,
                                   morphing, inlet, rho, t, dt);
      const Eigen::Vector3d n = vb.at(i, j).Normal();
      const double dS = vb.at(i, j).CalcB() * vb.at(i, j).CalcC();

      res += n * deltaP * dS;

      drag += CalcDragOnPanel(i, j, vb, vb_prev, wake_first, wake_last,
                              morphing, inlet, rho, t, dt);
    }
  }
  res.x() = drag;
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

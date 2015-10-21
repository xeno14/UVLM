/**
 * @file katz_plotkin.cpp
 * @brief Add description here
 */

#include "../util.h"
#include "katz_plotkin.h"


namespace UVLM {
namespace calc_load {
namespace internal {

}  // namespace internal

AerodynamicLoad CalcLoadKatzPlotkin(const VortexContainer& c, const VortexContainer& cp,
                         const Morphing& morphing, const UVLMVortexRing& rings,
                         const Eigen::Vector3d& freestream, const double rho,
                         const double t, const double dt) {
  
  LOG(INFO) << "Katz";
  Eigen::Vector3d F;
  double Pin=0, Pout=0;

  auto dim = DoubleLoop(c.rows(), c.cols());
  for (std::size_t index=0; index<dim.size(); ++index) {
    const auto i = dim[index].first, j = dim[index].second;
    const auto& v = c.at(dim[index].first, dim[index].second);

    const Eigen::Vector3d centroid = v.Centroid();
    const double db = v.CalcB();    // panel length
    const double dc = v.CalcC();    // panel length

    // prepare
    const auto Um = internal::CalcUm(morphing, centroid, freestream, t);
    const auto P = internal::CalcProjectionOperator(Um);
    const auto n = v.Normal();
    const double alpha = v.AngleOfAttack(Um);

    Eigen::Vector3d Uw, Ubc;
    rings.InducedVelocityByWake(&Uw, centroid);
    rings.InducedVelocityByChordwiseBound(&Ubc, centroid);

    // steady part
    const double Lst =
        internal::CalcLocalLiftSt(Um, Uw, c.Grad(i, j), db, dc, alpha, rho);
    const double Dst =
        internal::CalcLocalDragSt(Um, Ubc, P, n, c.DeltaGamma(i, j), db, rho);

    // unstready part
    double dg_dt = (c.at(i, j).gamma() - cp.at(i, j).gamma()) /dt;
    if (i > 0) {
      dg_dt += (c.at(i - 1, j).gamma() - cp.at(i - 1, j).gamma()) / dt;
    }
    dg_dt /= 2;
    const double Lunst = rho * dg_dt * db * dc * cos(alpha);
    const double Dunst = rho * dg_dt * db * dc * sin(alpha);

    // sum
    const Eigen::Vector3d e_drag = Um / Um.norm();
    const Eigen::Vector3d e_lift = P * n;
    F += e_drag * (Dst + Dunst) + e_lift * (Lst + Lunst);
  }

  return AerodynamicLoad{F, Pin, Pout};
}


}  // namespace calc_load
}  // namespace UVLM

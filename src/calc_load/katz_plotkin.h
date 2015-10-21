
/**
 * @file katz_plotkin.h
 * @brief Add description here
 */
#pragma once

#include "functions.h"

#include "../../proto/uvlm.pb.h"
#include "../vortex_container.h"
#include "../shed.h"
#include "../morphing.h"
#include "../uvlm_vortex_ring.h"
#include "functions.h"

#include <glog/logging.h>

namespace UVLM {
namespace calc_load {
namespace internal {
inline double CalcLocalLiftSt(const Eigen::Vector3d& Um, const Eigen::Vector3d& Uw,
                     const Eigen::Vector3d& grad, const double A,
                     const double alpha, const double rho) {
  return rho * (Um + Uw).dot(grad) * A * cos(alpha);
}

inline double CalcLocalDragSt(const Eigen::Vector3d& Ubc,
                              const Eigen::Vector3d& Uw,
                              const Eigen::Matrix3d& P,
                              const Eigen::Vector3d& n, const double dg,
                              const double b, const double rho) {
  return rho * (-(Ubc + Uw).dot(P * n) * dg * b);
}

// TODO unsteady term

}  // namespace internal

AerodynamicLoad CalcLoadKatzPlotkin(const VortexContainer& c,
                                    const VortexContainer& cp,
                                    const Morphing& morphing,
                                    const UVLMVortexRing& rings,
                                    const Eigen::Vector3d& freestream,
                                    const double rho, const double t,
                                    const double dt);
}}

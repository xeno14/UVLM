
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

inline double CalcLocalLift(const Eigen::Vector3d& Um, const Eigen::Vector3d& Uw,
                     const Eigen::Vector3d& grad, double dg_dt, const double A,
                     const double alpha, const double rho) {
  return rho * ((Um + Uw).dot(grad) + dg_dt) * A * cos(alpha);
}

}}

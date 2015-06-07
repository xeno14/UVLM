/**
 * @file proto_adaptor.h
 * @brief Add description here
 */
#pragma once

#include "../proto/uvlm.pb.h"
#include "uvlm_vortex_ring.h"


namespace UVLM {

inline Eigen::Vector3d PointToVector3d(const proto::Point& point) {
  return Eigen::Vector3d(point.x(), point.y(), point.z());
}

}  // namespace UVLM

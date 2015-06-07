/**
 * @file proto_adaptor.h
 * @brief Add description here
 */
#pragma once

#include "../proto/uvlm.pb.h"
#include "uvlm_vortex_ring.h"


namespace UVLM {

inline void PointToVector3d(Eigen::Vector3d* result,
                            const proto::Point& point) {
  *result << point.x() << point.y() << point.z();
}

}  // namespace UVLM

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

inline proto::Point Vector3dToPoint(const Eigen::Vector3d& vec) {
  proto::Point res;
  res.set_x(vec.x());
  res.set_y(vec.y());
  res.set_z(vec.z());
  return res;
}

inline proto::VortexRing VortexRingToProto(const VortexRing& v) {
  proto::VortexRing res;
  res.set_gamma(v.gamma());
  for (const auto& node : v.nodes()) {
    auto* target = res.add_nodes();
    *target = Vector3dToPoint(node);
  }
  for (const auto& node : v.nodes0()) {
    auto* target = res.add_nodes0();
    *target = Vector3dToPoint(node);
  }
  return res;
}

inline void UVLMVortexRingToBird(proto::FlyingWing* bird,
                                 const UVLMVortexRing& rings) {
  bird->Clear();
  for (const auto& v : rings.bound_vortices()) {
    auto* target = bird->add_bound_vortices();
    *target = VortexRingToProto(v);
  }
  for (const auto& v : rings.wake_vortices()) {
    auto* target = bird->add_wake_vortices();
    *target = VortexRingToProto(v);
  }
  auto* origin = bird->mutable_origin();
  *origin = Vector3dToPoint(rings.origin());
}

}  // namespace UVLM

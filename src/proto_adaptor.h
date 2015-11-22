/**
 * @file proto_adaptor.h
 * @brief Adaptors from proto to UVLM classes
 */
#pragma once

#include "../proto/uvlm.pb.h"
#include "vortex_container.h"
#include "wing_builder.h"
#include "uvlm_vortex_ring.h"

namespace UVLM {

/** @brief proto::Point to Eigen::Vector3d */
inline Eigen::Vector3d PointToVector3d(const proto::Point& point) {
  return Eigen::Vector3d(point.x(), point.y(), point.z());
}

/** @brief Eigen::Vector3d to proto::Point */
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

/**
 *
 * @todo node0のコピー
 */
inline UVLM::VortexRing ProtoToVortexRing(const proto::VortexRing& v) {
  UVLM::VortexRing res;
  res.set_gamma(v.gamma());
  for (const auto& node : v.nodes()) {
    res.PushNode(PointToVector3d(node));
  }
  return res;
}

template <class Range>
inline std::vector<Eigen::Vector3d> PointsToVector(const Range& points) {
  std::vector<Eigen::Vector3d> res;
  std::transform(points.begin(), points.end(), std::back_inserter(res),
                 [](const auto& p) { return UVLM::PointToVector3d(p); });
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

inline auto Snapshot2ToContainers(
    std::vector<UVLM::VortexContainer>* containers,
    const proto::Snapshot2& snapshot2) {
  auto vortices = std::make_shared<std::vector<UVLM::VortexRing>>(
      snapshot2.vortices().size());
  std::transform(snapshot2.vortices().begin(), snapshot2.vortices().end(),
                 vortices->begin(), ProtoToVortexRing);

  containers->clear();
  for (const auto& shape : snapshot2.container_shapes()) {
    containers->emplace_back(vortices, shape.rows(), shape.cols(), shape.id(),
                             shape.chord(), shape.span());
  }
  return vortices;
}

inline void Snapshot2ToMorphingVelocities(
    std::vector<Eigen::Vector3d>* centers,
    std::vector<std::vector<Eigen::Vector3d>>* v_nodes,
    std::vector<Eigen::Vector3d>* freestreams,
    const proto::Snapshot2& snapshot2) {
  centers->clear();
  v_nodes->clear();
  for (const auto& mv : snapshot2.morphing_velocities()) {
    centers->emplace_back(PointToVector3d(mv.center()));
    freestreams->emplace_back(PointToVector3d(mv.freestream()));
    std::vector<Eigen::Vector3d> a;
    for (const auto& vn : mv.nodes()) {
      a.emplace_back(PointToVector3d(vn));
    }
    v_nodes->emplace_back(a);
  }
}

}  // namespace UVLM

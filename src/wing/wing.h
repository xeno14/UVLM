
/**
 * @file wing_main.h
 * @brief Add description here
 */
#pragma once

#include "../../proto/uvlm.pb.h"
#include "../proto_adaptor.h"
#include "rect.h"
#include "naca00XX.h"
#include "naca4digit.h"

#include <eigen3/Eigen/Core>

namespace UVLM {
namespace wing {

inline void SetOrigin(::UVLM::proto::Wing* wing,
                      const Eigen::Vector3d& origin) {
  auto* op = wing->mutable_origin();
  op->set_x(origin.x());
  op->set_y(origin.y());
  op->set_z(origin.z());
}

void WholeWing(::UVLM::proto::Wing* wing, const ::UVLM::proto::Wing& half);

template <class... Args>
std::unique_ptr<WingGenerator> GeneratorFactory(const std::string& name,
                                                Args... args) {
  // TODO
  return nullptr;
}

}  // namespace wing
}  // namespace UVLM

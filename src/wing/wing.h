
/**
 * @file wing_main.h
 * @brief Add description here
 */
#pragma once

#include "rect.h"
#include "naca00XX.h"
#include "../../proto/uvlm.pb.h"
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

}  // namespace wing
}  // namespace UVLM

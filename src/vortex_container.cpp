/**
 * @file vortex_container.cpp
 * @brief Add description here
 */

#include "vortex_container.h"

#include <glog/logging.h>

namespace UVLM {

// TODO implement
std::vector<Eigen::Vector3d> VortexContainer::DumpPos() const {
  std::vector<Eigen::Vector3d> res((cols()+1) * (rows()+1));
  return res;
}

// TODO implement
void VortexContainer::LoadPos(const std::vector<Eigen::Vector3d>& pos) {
  CHECK(pos.size() == (cols()+1) * (rows()+1)); 
}

}  // namespace UVLM

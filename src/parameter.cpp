/**
 * @file parameter.cpp
 * @brief Add description here
 */

#include "parameter.h"

namespace UVLM {
namespace parameter {

void InitParam(const YAML::Node& node) {
  internal::ParamManager::instance().InitParam(node);
}

}  // namespace parameter
}  // namespace UVLM

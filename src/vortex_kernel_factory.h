
/**
 * @file vortex_kernel_factory.h
 * @brief Add description here
 */
#pragma once

#include "vortex_kernel.h"

#include <algorithm>
#include <boost/algorithm/string/join.hpp>
#include <glog/logging.h>
#include <map>

namespace {

enum VortexKernelType {
  CUTOFF,
  ROSENHEAD_MOORE,
};

const std::map<std::string, VortexKernelType> NAME_MAP{
    {"cutoff", VortexKernelType::CUTOFF},
    {"rosenhead-moore", VortexKernelType::ROSENHEAD_MOORE},
};

auto GetKeys() {
  std::vector<std::string> keys;
  std::transform(NAME_MAP.begin(), NAME_MAP.end(), std::back_inserter(keys),
                 [](const auto& e) { return e.first; });
  return keys;
}

const std::string NAMES = boost::algorithm::join(GetKeys(), ", ");

}  // anonymous namespace
namespace UVLM {
namespace vortex_kernel {

template <class... Args>
VortexKernel* VortexKernelFactory(const std::string& name, double a,
                                  Args... args) {
  if (NAME_MAP.find(name) == NAME_MAP.end()) {
    LOG(FATAL) << "Invalid advection name! Please choose from following: ["
               << NAMES << "]";
  }
  switch (NAME_MAP.at(name)) {
    case CUTOFF:
      return new CutOffKernel(a);
    case ROSENHEAD_MOORE:
      return new RosenheadMooreKernel(a);
    default:
      return nullptr;
  }
}

}  // namespace vortex_kernel
}  // namespace UVLM

/**
 * @file advect_factory.cpp
 * @brief Add description here
 */

#include "advect_factory.h"

#include <algorithm>
#include <boost/algorithm/string/join.hpp>
#include <glog/logging.h>
#include <map>

namespace {

enum AdvectionType {
  EULER,
  RK2,
  RK4,
  AB2,
};

const std::map<std::string, AdvectionType> NAME_MAP{
    {"euler", AdvectionType::EULER},
    {"RK2", AdvectionType::RK2},
    {"RK4", AdvectionType::RK4},
    {"AB2", AdvectionType::AB2},
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
namespace advect {

Advection* AdvectFactory(const std::string& name) {
  if (NAME_MAP.find(name) == NAME_MAP.end()) {
    LOG(FATAL) << "Invalid advection name! Please choose from following: ["
               << NAMES << "]";
  }
  switch(NAME_MAP.at(name)) {
    case EULER:
      return new Euler;
    case RK2:
      return new RungeKutta2;
    case RK4:
      return new RungeKutta4;
    case AB2:
      return new AdamsBashforth2;
    default:
      return nullptr;
  }
}

std::string GetNames() { return NAMES; }

}  // namespace advect
}  // namespace UVLM

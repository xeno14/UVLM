/**
 * @file calc_load.cpp
 * @brief Add description here
 */

#include "calc_load.h"

namespace UVLM {

namespace calc_load {

double CalcC (const VortexRing& v) {
  return (v.nodes()[0] - v.nodes()[1]).norm();
}

double CalcB (const VortexRing& v) {
  return (v.nodes()[0] - v.nodes()[3]).norm();
}

}  // namespace internal
}  // namespace UVLM

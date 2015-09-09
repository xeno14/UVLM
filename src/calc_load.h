
/**
 * @file calc_load.h
 * @brief Add description here
 */
#pragma once

#include "vortex_container.h"

namespace UVLM {

namespace calc_load {

/**
 * 渦輪のchord方向の長さを求める
 */
double CalcC(const VortexRing& v);

/**
 * 渦輪のspan方向の長さを求める
 */
double CalcB(const VortexRing& v);

}  // namespace internal

}  // namespace UVLM

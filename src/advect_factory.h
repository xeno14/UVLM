
/**
 * @file advect_factory.h
 * @brief Add description here
 */
#pragma once

#include "advect.h"

namespace UVLM {
namespace advect {

Advection* AdvectFactory(const std::string& name);

std::string GetNames();

}  // namespace advect
}  // namespace UVLM

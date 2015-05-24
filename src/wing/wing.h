
/**
 * @file wing_main.h
 * @brief Add description here
 */
#pragma once

#include <iostream>

namespace UVLM {
namespace wing {

inline std::ostream& output(std::ostream& os, double x, double y, double z) {
  return os << x << "\t" << y << "\t" << z << std::endl;
}

void Perfome(std::ostream& os);

}
}

/**
 * @file naca00XX.cpp
 * @brief Add description here
 */

#include "naca00XX.h"

#include <gflags/gflags.h>
#include <cmath>
#include <iostream>

namespace UVLM {
namespace wing {

double NACA00XX(double x, double c, int xx) {
  double t = xx * 0.01;
  double z = x / c;
  return 5 * t * c * (0.2969 * sqrt(z) + (-0.1260) * (z) + (-0.3516) * z * z +
                      0.2843 * z * z * z + (-0.1015) * z * z * z * z);
}

void NACA00XXGenerator::Generate(UVLM::proto::Wing* wing, const double chord,
                                 const double span, const std::size_t rows,
                                 const std::size_t cols) const {
  // 長方形の翼を生成してzの値のみを変更する
  RectGenerator::Generate(wing, chord, span, rows, cols);

  for (int i = 0; i < wing->points_size(); i++) {
    auto* p = wing->mutable_points(i);
    p->set_z(NACA00XX(p->x(), chord, digit_));
  }
}

}  // namespace wing
}  // namespace UVLM

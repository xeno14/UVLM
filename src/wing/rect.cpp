/**
 * @file rect.cpp
 * @brief Add description here
 */

#include "rect.h"

namespace UVLM {
namespace wing {

void RectGenerator::Generate(UVLM::proto::Wing* wing) {
  WingGenerator::Generate(wing);

  const double chord = chord_;
  const double span = span_;
  const int rows = rows_;
  const int cols = cols_;
  const double dx = chord / rows;
  const double dy = span / cols;

  for (int i = 0; i <= rows; i++) {
    for (int j = 0; j <= cols; j++) {
      double x = i * dx;
      double y = j * dy;
      double z = 0;
      auto* point = wing->add_points();
      point->set_x(x);
      point->set_y(y);
      point->set_z(z);
    }
  }
}

}  // namespace wing
}  // namespace UVLM

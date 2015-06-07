/**
 * @file naca00XX.cpp
 * @brief Add description here
 */

#include "naca00XX.h"
#include "wing.h"

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

void NACA00XXGenerator::Generate(UVLM::proto::Wing* wing) {
  const int xx = digit_;
  const double chord = chord_;
  const double span = span_;
  const int rows = rows_;
  const int cols = cols_;
  const double dx = chord / rows;
  const double dy = span / cols;

  if (verbose_) {
    std::cerr << "NACA00" << xx << std::endl;
    std::cerr << "Chord: " << chord << std::endl;
    std::cerr << "Span: " << span << std::endl;
    std::cerr << rows << "x" << cols << std::endl;
    std::cerr << dx << "@" << dy << "\n";
  }

  wing->set_cols(cols);
  wing->set_rows(rows);
  for (int i = 0; i <= rows; i++) {
    for (int j = 0; j <= cols; j++) {
      double x = i * dx;
      double y = j * dy;
      double z = NACA00XX(x, chord, xx);
      auto* point = wing->add_points();
      point->set_x(x);
      point->set_y(y);
      point->set_z(z);
    }
  }
}

}  // namespace wing
}  // namespace UVLM

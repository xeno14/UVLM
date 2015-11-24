/**
 * @file naca4digit.cpp
 * @brief Add description here
 */

#include "naca4digit.h"


namespace UVLM {
namespace wing {

double NACA4digit(double x, double c, int xx) {
  int first = xx / 10;
  int second = xx % 10;
  const double m = first / 100.0;
  const double p = second / 10.0;

  double res = 0;
  if (0 <= x && x <= p * c) {
    res = m * x / p / p * (2 * p - x / c);
  } else if (p * c <= x && x <= c) {
    res = m * (c - x) / (1 - p) / (1 - p) * (1 + x / c - 2 * p);
  }
  return res;
}

void NACA4digitGenerator::Generate(UVLM::proto::Wing* wing, const double chord,
                                   const double span, const std::size_t rows,
                                   const std::size_t cols) {
  // map cordinates on a rectangular wing
  RectGenerator::Generate(wing, chord, span, rows, cols);

  for (int i = 0; i < wing->points_size(); i++) {
    auto* p = wing->mutable_points(i);
    p->set_z(NACA4digit(p->x(), chord, digit_));
  }
}

}  // namespace wing
}  // namespace UVLM

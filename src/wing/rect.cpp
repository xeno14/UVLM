/**
 * @file rect.cpp
 * @brief Add description here
 */

#include "rect.h"

#include <iostream>
#include <algorithm>
#include <gflags/gflags.h>

DEFINE_bool(chebyshev_chord, false, "use chebyshev in chordwise direction");
DEFINE_bool(chebyshev_span, false, "use chebyshev in spanwise direction");

namespace UVLM {
namespace wing {

void RectGenerator::Generate(UVLM::proto::Wing* wing, const double chord,
                        const double span, const std::size_t rows,
                        const std::size_t cols) const {
  WingGenerator::Generate(wing, chord, span, rows, cols);

  std::vector<double> xs, ys;
  if (FLAGS_chebyshev_chord) {
    TransformChebyshev(&xs, rows + 1, 0, chord);
  } else {
    xs = linspace(0, chord, rows + 1);
  }
  if (FLAGS_chebyshev_span) {
    TransformChebyshev(&ys, cols + 1, 0, span);
  } else {
    ys = linspace(0, span, cols + 1);
  }

  for (const double x : xs) {
    for (const double y : ys) {
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

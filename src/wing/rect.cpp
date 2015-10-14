/**
 * @file rect.cpp
 * @brief Add description here
 */

#include "../util.h"
#include "rect.h"

#include <algorithm>
#include <gflags/gflags.h>

DEFINE_bool(chebyshev_chord, false, "use chebyshev in chordwise direction");
DEFINE_bool(chebyshev_span, false, "use chebyshev in spanwise direction");

namespace UVLM {
namespace wing {

void RectGenerator::Generate(UVLM::proto::Wing* wing) {
  WingGenerator::Generate(wing);

  const double chord = chord_;
  const double span = span_;
  const int rows = rows_;
  const int cols = cols_;

  std::vector<double> xs(rows), ys(cols);
  if (FLAGS_chebyshev_chord) {
    auto theta = linspace(M_PI_2, 0, rows + 1);
    std::transform(theta.begin(), theta.end(), xs.begin(),
                   [chord](double t) { return chord * cos(t); });
  } else {
    xs = linspace(0, chord, rows + 1);
  }
  if (FLAGS_chebyshev_span) {
    auto theta = linspace(M_PI_2, 0, cols + 1);
    std::transform(theta.begin(), theta.end(), ys.begin(),
                   [span](double t) { return span * cos(t); });
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

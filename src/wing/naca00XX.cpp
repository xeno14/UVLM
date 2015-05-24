/**
 * @file naca00XX.cpp
 * @brief Add description here
 */

#include "wing.h"

#include <gflags/gflags.h>

#include <cmath>

DEFINE_int32(XX, 12, "Digit of NACA00XX");
DEFINE_double(chord, 1.0, "Chord length");
DEFINE_double(aspect_ratio, 4.0, "Aspect ratio of the wing");
DEFINE_int32(rows, 10, "Number of rows (x).");


namespace UVLM {
namespace wing {

/**
 * @brief Equation for a symmetrical 4-digit NACA airfoil
 * http://en.wikipedia.org/wiki/NACA_airfoil#Equation_for_a_symmetrical_4-digit_NACA_airfoil
 * @param c chord length
 * @param xx number
 */
inline double NACA00XX(double x, double c, int xx) {
  double t = xx * 0.01;
  double z = x / c;
  return 5 * t * c * (0.2969 * sqrt(z) + (-0.1260) * (z) + (-0.3516) * z * z +
                      0.2843 * z * z * z + (-0.1015) * z * z * z * z);
}

void Perfome(std::ostream& os) {
  const int xx = FLAGS_XX;
  const double chord = FLAGS_chord;
  const double span = chord * FLAGS_aspect_ratio;
  const int rows = FLAGS_rows;
  const double dx = chord/rows;
  const int cols = int(span/dx);
  const double dy = span / cols;
  std::cerr << "NACA00" << xx << std::endl;
  std::cerr << "Chord: " << chord << std::endl;
  std::cerr << "Span: " << span<< std::endl;
  std::cerr << rows << "x" << cols << std::endl;
  std::cerr << dx << "@" << dy <<"\n";
  for (int i=0; i<=rows; i++) {
    for (int j=0; j<=cols; j++) {
      double x = i * dx;
      double y = -span/2 + j * dy;
      double z = NACA00XX(x, chord, xx);
      output(os, x, y, z);
    }
  }
}

}  // namespace wing
}  // namespace UVLM

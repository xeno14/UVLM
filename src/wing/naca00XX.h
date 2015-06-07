
/**
 * @file naca00XX.h
 * @brief Add description here
 */
#pragma once

#include "wing_generator.h"

namespace UVLM {
namespace wing {

/**
 * @brief Equation for a symmetrical 4-digit NACA airfoil
 * http://en.wikipedia.org/wiki/NACA_airfoil#Equation_for_a_symmetrical_4-digit_NACA_airfoil
 * @param c chord length
 * @param xx number
 */
double NACA00XX(double x, double c, int xx);

class NACA00XXGenerator : public WingGenerator {
 public:
  NACA00XXGenerator(int digit, double chord, double ar, std::size_t rows,
                    std::size_t cols)
      : WingGenerator(rows, cols),
        digit_(digit),
        chord_(chord),
        aspect_ratio_(ar) {}
  void Generate(UVLM::proto::Wing* wing);

 private:
  int digit_;
  double chord_, aspect_ratio_;
};

}  // namespace wing
}  // namespace UVLM

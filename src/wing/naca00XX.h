
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
  NACA00XXGenerator(int digit, double chord, double span, std::size_t rows,
                    std::size_t cols)
      : WingGenerator(rows, cols),
        digit_(digit),
        chord_(chord),
        span_(span),
        verbose_(false) {}
  void Generate(UVLM::proto::Wing* wing);
  void set_verbose(bool flag) { verbose_ = flag; }

 private:
  int digit_;
  double chord_, span_;
  bool verbose_;
};

}  // namespace wing
}  // namespace UVLM

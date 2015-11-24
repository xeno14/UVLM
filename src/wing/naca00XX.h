
/**
 * @file naca00XX.h
 * @brief Add description here
 */
#pragma once

#include "rect.h"

namespace UVLM {
namespace wing {

/**
 * @brief Equation for a symmetrical 4-digit NACA airfoil
 * http://en.wikipedia.org/wiki/NACA_airfoil#Equation_for_a_symmetrical_4-digit_NACA_airfoil
 * @param c chord length
 * @param xx number
 */
double NACA00XX(double x, double c, int xx);

class NACA00XXGenerator : public RectGenerator {
 public:
  NACA00XXGenerator(int digit)
      : RectGenerator(),
        digit_(digit),
        verbose_(false) {}
  virtual ~NACA00XXGenerator() = default;
  virtual void Generate(UVLM::proto::Wing* wing,
                        const double chord, const double span,
                        const std::size_t rows,
                        const std::size_t cols) const override;
  void set_verbose(bool flag) { verbose_ = flag; }

 private:
  int digit_;
  bool verbose_;
};

}  // namespace wing
}  // namespace UVLM

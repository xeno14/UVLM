/**
 * @file naca4digit.h
 * @brief Add description here
 */
#pragma once


#include "rect.h"

namespace UVLM {
namespace wing {

/**
 * @brief mean camber of 4-digit NACA airfoil
 * Because the thickness is irrelavant, only first two digits are used.
 *
 * @param c chord length
 * @param xx number
 */
double NACA4digit(double x, double c, int xx);

class NACA4digitGenerator : public RectGenerator {
 public:
  NACA4digitGenerator(int digit, double chord, double span, std::size_t rows,
                    std::size_t cols)
      : RectGenerator(chord, span, rows, cols),
        digit_(digit),
        verbose_(false) {}
  virtual ~NACA4digitGenerator() = default;
  virtual void Generate(UVLM::proto::Wing* wing) override;
  void set_verbose(bool flag) { verbose_ = flag; }

 private:
  int digit_;
  bool verbose_;
};

}  // namespace wing
}  // namespace UVLM


/**
 * @file rect.h
 * @brief Add description here
 */
#pragma once

#include "../util.h"
#include "wing_generator.h"

#include <algorithm>
#include <vector>

namespace UVLM {
namespace wing {

inline void TransformChebyshev(std::vector<double>* result, std::size_t N,
                               double xmin, double xmax) {
  auto theta = linspace(-M_PI_2, 0, N);
  result->clear();
  std::transform(
      theta.cbegin(), theta.cend(), std::back_inserter(*result),
      [xmin, xmax](double t) { return xmin + (xmax - xmin) * cos(t); });
  *result->begin() = xmin;
  *result->rbegin() = xmax;
}

class RectGenerator : public WingGenerator {
 public:
  RectGenerator() : WingGenerator() {}
  virtual ~RectGenerator() = default;
  virtual void Generate(UVLM::proto::Wing* wing, const double chord,
                        const double span, const std::size_t rows,
                        const std::size_t cols) override;
};

}  // namespace wing
}  // namespace UVLM

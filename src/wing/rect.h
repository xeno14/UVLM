
/**
 * @file rect.h
 * @brief Add description here
 */
#pragma once

#include "wing_generator.h"

namespace UVLM {
namespace wing {

class RectGenerator : public WingGenerator {
 public:
  RectGenerator(double chord, double span, std::size_t rows, std::size_t cols)
      : WingGenerator(chord, span, rows, cols) {}
  virtual ~RectGenerator() = default;
  virtual void Generate(UVLM::proto::Wing* wing) override;
};

}  // namespace wing
}  // namespace UVLM

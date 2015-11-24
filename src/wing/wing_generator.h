
/**
 * @file wing_generator.h
 * @brief Add description here
 */
#pragma once

#include "../../proto/uvlm.pb.h"

namespace UVLM {
namespace wing {

class WingGenerator {
 public:
  WingGenerator() {}
  virtual ~WingGenerator() = default;

  virtual void Generate(UVLM::proto::Wing* wing, const double chord,
                        const double span, const std::size_t rows,
                        const std::size_t cols) {
    wing->Clear();
    wing->set_chord(chord);
    wing->set_span(span);
    wing->set_cols(cols);
    wing->set_rows(rows);
  }

  virtual void operator()(UVLM::proto::Wing* wing, const double chord,
                          const double span, const std::size_t rows,
                          const std::size_t cols) final {
    this->Generate(wing, chord, span, rows, cols);
  }
};

}  // namespace wing
}  // namespace UVLM

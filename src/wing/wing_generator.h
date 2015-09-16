
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
  WingGenerator(double chord, double span, std::size_t rows, std::size_t cols)
      : chord_(chord), span_(span), rows_(rows), cols_(cols) {}
  virtual ~WingGenerator() = default;
  virtual void Generate(UVLM::proto::Wing* wing) {
    wing->Clear();
    wing->set_chord(chord_);
    wing->set_span(span_);
    wing->set_cols(cols_);
    wing->set_rows(rows_);
  }

  virtual void operator() (UVLM::proto::Wing* wing) final {
    this->Generate(wing);
  }
 protected:
  double chord_, span_;
  std::size_t rows_, cols_;
};

}  // namespace wing
}  // namespace UVLM

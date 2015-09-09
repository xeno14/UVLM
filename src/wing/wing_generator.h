
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
  WingGenerator(std::size_t rows, std::size_t cols)
      : rows_(rows), cols_(cols) {}
  virtual ~WingGenerator() = default;
  virtual void Generate(UVLM::proto::Wing* wing) = 0;
  virtual void operator() (UVLM::proto::Wing* wing) final {
    this->Generate(wing);
  }
 protected:
  std::size_t rows_, cols_;
};

}  // namespace wing
}  // namespace UVLM

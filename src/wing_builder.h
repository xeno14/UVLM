/**
 * @file wing_builder.h
 * @brief Add description here
 */
#pragma once

#include "../proto/uvlm.pb.h"
#include "vortex_container.h"

namespace UVLM {
namespace internal {

struct WingHolder {
  std::vector<Eigen::Vector3d> points;
  Eigen::Vector3d origin;
  std::size_t rows, cols;
};

template <class InputIterator>
std::vector<Eigen::Vector3d> TransfromRow(InputIterator first,
                                          InputIterator last) {
  std::vector<Eigen::Vector3d> row(first, last);
  std::vector<Eigen::Vector3d> reversed_row = row;

  std::reverse(reversed_row.begin(), reversed_row.end());
  std::for_each(reversed_row.begin(), reversed_row.end(),
                [](Eigen::Vector3d& pos) { pos.y() *= -1; });

  std::vector<Eigen::Vector3d> res;
  res.insert(res.end(), reversed_row.begin(), reversed_row.end());
  res.insert(res.end(), row.begin() + 1, row.end());
  return res;
}

}  // namespace internal

class WingBuilder {
 public:
  WingBuilder(std::vector<VortexContainer>* containers,
              std::shared_ptr<std::vector<VortexRing>> vortices)
      : containers_(containers), vortices_(vortices), is_built_(false) {}

  ~WingBuilder() {
    if (!is_built_) Build();
  }

  WingBuilder& AddWing(const proto::Wing& wing);

  void Build();

  void BuildWing(VortexContainer* container,
                 const internal::WingHolder holder) const;

  static std::size_t CountTotalSize(
      const std::vector<internal::WingHolder>& holders_);

  const std::vector<internal::WingHolder>& holders() const { return holders_; }

 private:
  std::vector<VortexContainer>* containers_;
  std::shared_ptr<std::vector<VortexRing>> vortices_;
  std::vector<internal::WingHolder> holders_;
  bool is_built_;
};

}  // namespace UVLM

/**
 * @file wing_builder.cpp
 * @brief Add description here
 */

#include "wing_builder.h"
#include "proto_adaptor.h"

namespace UVLM {

WingBuilder& WingBuilder::AddWing(const proto::Wing& wing) {
  const std::size_t point_rows = wing.rows() + 1;
  const std::size_t point_cols = wing.cols() + 1;

  internal::WingHolder holder;
  holder.rows = wing.rows();
  holder.cols = wing.cols() * 2;

  // ランダムアクセス可能なvectorに変形する
  std::vector<Eigen::Vector3d> points(wing.points().size());
  std::transform(wing.points().begin(), wing.points().end(), points.begin(),
                 PointToVector3d);
  for (std::size_t i=0; i < point_rows; i++) {
    auto first = points.begin() + i * point_cols;
    auto last = first + point_cols;
    auto row = internal::TransfromRow(first, last);
    holder.points.insert(holder.points.end(), row.begin(), row.end());
  }
  holders_.push_back(holder);
  return *this;
}

void WingBuilder::Build() {
  vortices_->resize(CountTotalSize(holders_));
}

std::size_t WingBuilder::CountTotalSize(
      const std::vector<internal::WingHolder>& holders_) {
  std::size_t res = 0;
  for (const auto& holder : holders_) res += holder.rows * holder.cols;
  return res;
}

}  // namespace UVLM

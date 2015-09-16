/**
 * @file wing_builder.cpp
 * @brief Add description here
 */

#include "wing_builder.h"
#include "proto_adaptor.h"

#include <iostream>

namespace UVLM {

WingBuilder& WingBuilder::AddWing(const proto::Wing& wing) {
  is_modified_ = true;

  const std::size_t point_rows = wing.rows() + 1;
  const std::size_t point_cols = wing.cols() + 1;

  internal::WingHolder holder;
  holder.rows = wing.rows();
  holder.cols  = wing.cols() * 2;
  holder.chord = wing.chord();
  holder.span  = wing.span();
  if (wing.has_origin()) holder.origin = PointToVector3d(wing.origin());
  else holder.origin = Eigen::Vector3d::Zero();

  // ランダムアクセス可能なvectorに変形する
  std::vector<Eigen::Vector3d> points;
  std::transform(wing.points().begin(), wing.points().end(),
                 std::back_inserter(points), PointToVector3d);

  for (std::size_t i = 0; i < point_rows; i++) {
    auto first = points.begin() + i * point_cols;
    auto last = first + point_cols;
    auto row = internal::TransformRow(first, last);
    holder.points.insert(holder.points.end(), row.begin(), row.end());
  }

  // 原点をずらす
  for (auto& point : holder.points) point += holder.origin;

  holders_.push_back(holder);
  return *this;
}

void WingBuilder::Build() {
  if (!is_modified_) return;

  vortices_->resize(CountTotalSize(holders_));
  containers_->resize(holders_.size());

  for (std::size_t i = 0; i < containers_->size(); i++) {
    containers_->at(i).set_vortices(vortices_, holders_[i].rows,
                                    holders_[i].cols, i, holders_[i].chord,
                                    holders_[i].span);
    BuildWing(&containers_->at(i), holders_[i]);
  }

  is_modified_ = false;
}

void WingBuilder::BuildWing(VortexContainer* container,
                            const internal::WingHolder& holder) {
  const std::size_t point_cols = holder.cols + 1;
  for (std::size_t i = 0; i < holder.rows; i++) {
    for (std::size_t j = 0; j < holder.cols; j++) {
      std::size_t indices[4] = {
          j + i * point_cols,            // 0
          j + (i + 1) * point_cols,      // 1
          j + 1 + (i + 1) * point_cols,  // 2
          j + 1 + i * point_cols,        // 3
      };
      VortexRing vortex;
      for (auto idx : indices) vortex.nodes().push_back(holder.points[idx]);
      vortex.SaveReferenceNode();
      container->at(i, j) = vortex;
    }
  }
}

std::size_t WingBuilder::CountTotalSize(
      const std::vector<internal::WingHolder>& holders_) {
  std::size_t res = 0;
  for (const auto& holder : holders_) res += holder.rows * holder.cols;
  return res;
}

}  // namespace UVLM

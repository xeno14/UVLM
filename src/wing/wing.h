
/**
 * @file wing_main.h
 * @brief Add description here
 */
#pragma once

#include "../../proto/uvlm.pb.h"
#include "../proto_adaptor.h"
#include "rect.h"
#include "naca00XX.h"
#include "naca4digit.h"

#include <eigen3/Eigen/Core>

namespace UVLM {
namespace wing {

inline void SetOrigin(::UVLM::proto::Wing* wing,
                      const Eigen::Vector3d& origin) {
  auto* op = wing->mutable_origin();
  op->set_x(origin.x());
  op->set_y(origin.y());
  op->set_z(origin.z());
}

inline void WholeWing(::UVLM::proto::Wing* wing,
                      const ::UVLM::proto::Wing& half) {
  wing->Clear();
  wing->set_cols(half.cols() * 2);
  wing->set_rows(half.rows());
  wing->mutable_origin()->CopyFrom(half.origin());
  wing->set_chord(half.chord());
  wing->set_span(half.span());

  std::vector<::UVLM::proto::Point> half_points(half.points().begin(),
                                                half.points().end());
  std::vector<::UVLM::proto::Point> whole_points;
  std::vector<::UVLM::proto::Point> mirror;
  const std::size_t COLS = half.cols() + 1;
  const std::size_t ROWS = half.points().size() / COLS;
  for (std::size_t i=0; i<ROWS; ++i) {
    std::vector<::UVLM::proto::Point> row_points(
        half_points.begin() + i * COLS, half_points.begin() + (i + 1) * COLS);
    std::vector<::UVLM::proto::Point> mirror_row_points(COLS);
    std::transform(row_points.begin(), row_points.end(),
                   mirror_row_points.begin(),
                   [](const ::UVLM::proto::Point& p) {
                     ::UVLM::proto::Point res(p);
                     res.set_y(p.y() * (-1));
                     return res;
                   });
    whole_points.insert(whole_points.end(), mirror_row_points.rbegin(),
                        mirror_row_points.rend());
    whole_points.insert(whole_points.end(), row_points.begin()+1,
                        row_points.end());
  }
  for (const auto& p : whole_points) {
    wing->add_points()->CopyFrom(p);
  }
}

}  // namespace wing
}  // namespace UVLM

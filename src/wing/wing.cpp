/**
 * @file wing.cpp
 * @brief Add description here
 */

#include "wing.h"


namespace UVLM {
namespace wing {

void WholeWing(::UVLM::proto::Wing* wing, const ::UVLM::proto::Wing& half) {
  wing->Clear();
  wing->set_cols(half.cols() * 2);
  wing->set_rows(half.rows());
  wing->mutable_origin()->CopyFrom(half.origin());
  wing->set_chord(half.chord());
  wing->set_span(half.span());

  std::vector<Eigen::Vector3d> half_points =
      UVLM::PointsToVector(half.points());
  std::vector<Eigen::Vector3d> whole_points;
  std::vector<Eigen::Vector3d> mirror;

  const std::size_t COLS = half.cols() + 1;
  const std::size_t ROWS = half.points().size() / COLS;
  const Eigen::Vector3d origin = PointToVector3d(half.origin());

  for (std::size_t i = 0; i < ROWS; ++i) {
    std::vector<Eigen::Vector3d> row_points(
        half_points.begin() + i * COLS, half_points.begin() + (i + 1) * COLS);
    std::vector<Eigen::Vector3d> mirror_row_points(COLS);
    std::transform(row_points.begin(), row_points.end(),
                   mirror_row_points.begin(), [](const Eigen::Vector3d& p) {
                     Eigen::Vector3d res(p);
                     res.y() *= -1;
                     return res;
                   });
    whole_points.insert(whole_points.end(), mirror_row_points.rbegin(),
                        mirror_row_points.rend());
    whole_points.insert(whole_points.end(), row_points.begin() + 1,
                        row_points.end());
  }
  // move back origin
  for (auto& p : whole_points) p += origin;

  for (const auto& p : whole_points) {
    wing->add_points()->CopyFrom(Vector3dToPoint(p));
  }
}

}  // namespace wing
}  // namespace UVLM

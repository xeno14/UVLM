/**
 * @file simple_simulator.cpp
 * @brief Add description here
 */

#include "simple_simulator.h"
#include "../wing/wing.h"
#include "../output.h"

#include <glog/logging.h>

namespace UVLM {
namespace simulator {

std::vector<Eigen::Vector3d> CollocationPoints(
    const MultipleSheet<Eigen::Vector3d>& pos) {
  std::vector<Eigen::Vector3d> res;
  for (std::size_t n = 0; n < pos.num(); n++) {
    for (std::size_t i = 0; i < pos.rows() - 1; i++) {
      for (std::size_t j = 0; j < pos.cols() - 1; j++) {
        res.push_back((pos.at(n, i, j) + pos.at(n, i + 1, j) +
                       pos.at(n, i + 1, j + 1) + pos.at(n, i, j + 1)) /
                      4);
      }
    }
  }
  return res;
}

std::vector<Eigen::Vector3d> Normals(
    const MultipleSheet<Eigen::Vector3d>& pos) {
  std::vector<Eigen::Vector3d> res;
  for (std::size_t n = 0; n < pos.num(); n++) {
    for (std::size_t i = 0; i < pos.rows() - 1; ++i) {
      for (std::size_t j = 0; j < pos.cols() - 1; ++j) {
        Eigen::Vector3d nrml =
            (pos.at(n, i + 1, j + 1) - pos.at(n, i, j))
                .cross(pos.at(n, i, j + 1) - pos.at(n, i + 1, j));
        nrml.normalize();
        res.emplace_back(nrml);
      }
    }
  }
  return res;
}

Eigen::MatrixXd SimpleSimulator::CalcMatrix(const std::vector<Eigen::Vector3d>& cpos,
                           const std::vector<Eigen::Vector3d>& normal) const {
  // A_kl
  Eigen::MatrixXd res(wing_gamma_.size(), wing_gamma_.size());
  MultipleSheet<double> gamma(wing_gamma_);
  std::fill(gamma.begin(), gamma.end(), 1);

  // loop for all bound vortices
  for (auto index_K : wing_gamma_.list_index()) {
    std::size_t K;
    std::tie(K, std::ignore, std::ignore, std::ignore) = index_K;
    const auto& cp = cpos[K];
    const auto& nl = normal[K];

    for (auto index_L : wing_gamma_.list_index()) {
      std::size_t L, nn, ii, jj;
      std::tie(L, nn, ii, jj) = index_L;
      auto u = UVLM::VORING(cp, wing_pos_, gamma, nn, ii, jj);
      res(K, L) = u.dot(nl);
    }
  }
  return res;
}

void SimpleSimulator::AddWing(const Morphing& morphing, const double chord,
                              const double span, const std::size_t rows,
                              const std::size_t cols,
                              const Eigen::Vector3d& origin) {
  if (wing_info_.size()) {
    CHECK(wing_info_.rbegin()->rows == rows);
    CHECK(wing_info_.rbegin()->cols == cols);
  }
  wing_info_.push_back(
      WingInformation{morphing, chord, span, rows, cols, origin});
}

void SimpleSimulator::BuildWing() {
  CHECK(wing_info_.size() > 0) << "call AddWing at least once";
  const std::size_t num = wing_info_.size();
  const std::size_t rows = wing_info_.rbegin()->rows;
  const std::size_t cols = wing_info_.rbegin()->cols;

  wing_pos_.resize(num, rows + 1, cols + 1);
  wing_pos_init_.resize(num, rows + 1, cols + 1);
  wing_gamma_.resize(num, rows, cols);
  wing_gamma_prev_.resize(num, rows, cols);

  std::size_t n = 0;
  for (const auto& info : wing_info_) {
    const auto origin = info.origin;
    morphings_.push_back(info.morphing);
    morphings_.rbegin()->set_origin(origin);

    UVLM::proto::Wing wing, half;
    UVLM::wing::NACA4digitGenerator wing_generator(
        83, info.chord, info.span / 2, rows, cols / 2);
    wing_generator.Generate(&half);
    UVLM::wing::WholeWing(&wing, half);
    LOG(INFO) << "origin " << info.origin.transpose();
    auto points = UVLM::PointsToVector(wing.points());
    std::transform(points.begin(), points.end(),
        points.begin(),
        [origin](const auto& x) { return x + origin; });
    std::copy(points.begin(), points.end(), wing_pos_init_.sheet_begin(n));
    LOG(INFO) << wing_pos_.sheet_begin(n)->transpose();
    ++n;
  }
}

void SimpleSimulator::MainLoop(const std::size_t step, const double dt) {
  const double t = step * dt;

  std::copy(wing_gamma_.begin(), wing_gamma_.end(), wing_gamma_prev_.begin());
  
  LOG(INFO) << "Morphing";
  CHECK(wing_pos_.size() == wing_pos_init_.size());
  for (std::size_t n = 0; n < wing_pos_.num(); n++) {
    std::transform(wing_pos_init_.sheet_begin(n), wing_pos_init_.sheet_end(n),
                   wing_pos_.sheet_begin(n),
                   [t, n, this](const Eigen::Vector3d& x0) {
                     return this->morphings_[n].Perfome(x0, t);
                   });
  }
  
  const auto cpos = CollocationPoints(wing_pos_);
  const auto normal = Normals(wing_pos_);

  if (result_path_.size()) OutputPanels(step, dt);
}

void SimpleSimulator::Run(const std::size_t steps, const double dt) {
  BuildWing();
  for (std::size_t step = 1; step<=steps; step++) {
    MainLoop(step, dt);
  }
}

void SimpleSimulator::OutputPanels(const std::size_t step, const double dt) const {
  char filename[256];
  UVLM::proto::Snapshot2 snapshot;
  // LOG(INFO) << wing_gamma_.num();
  // LOG(INFO) << wing_pos_.num();
  for (std::size_t n = 0; n < wing_gamma_.num(); n++) {
    UVLM::output::SimpleAppendSnapshot(
        &snapshot, wing_pos_.sheet_begin(n), wing_pos_.sheet_end(n),
        wing_gamma_.sheet_begin(n), wing_gamma_.sheet_end(n),
        wing_gamma_.cols());
  }
  if (wake_gamma_.size()) {
    UVLM::output::SimpleAppendSnapshot(&snapshot, wake_pos_.begin(),
                                        wake_pos_.end(), wake_gamma_.begin(),
                                        wake_gamma_.end(), wake_gamma_.cols());
  }

  // TODO change path
  sprintf(filename, "%s/%08lu", result_path_.c_str(), step);
  std::ofstream ofs(filename);
  CHECK(ofs);
  snapshot.SerializeToOstream(&ofs);
}

}  // namespace simulator
}  // namespace UVLM

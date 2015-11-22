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

Eigen::Vector3d SimpleSimulator::BoundVelocity(const Eigen::Vector3d& x) const{
  Eigen::Vector3d res = Eigen::Vector3d::Zero();
  for (auto index : wing_gamma_.list_index()) {
    std::size_t n, i, j;
    std::tie(std::ignore, n, i, j) = index;
    res += UVLM::VORING(x, wing_pos_, wing_gamma_, n, i, j);
  }
  return res;
}

Eigen::Vector3d SimpleSimulator::WakeVelocity(const Eigen::Vector3d& x) const {
  Eigen::Vector3d res = Eigen::Vector3d::Zero();
  for (auto index : wake_gamma_.list_index()) {
    std::size_t n, i, j;
    std::tie(std::ignore, n, i, j) = index;
    res += UVLM::VORING(x, wake_pos_, wake_gamma_, n, i, j);
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

// auto SimpleSimulator::CalcRhs(const double t) const {
//   const std::size_t sz = cpos.size();
//   Eigen::VectorXd res(sz);
//   for (std::size_t K = 0; K < sz; ++K) {
//     Eigen::Vector3d u = -forward_flight + WakeVelocity(cpos[K]) -
//                         MorphingVelocity(cpos_init[K], t);
//     res(K) = -u.dot(normal[K]);
//   }
//   return res;
// }

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
    auto points = UVLM::PointsToVector(wing.points());
    std::transform(points.begin(), points.end(),
        points.begin(),
        [origin](const auto& x) { return x + origin; });
    std::copy(points.begin(), points.end(), wing_pos_init_.sheet_begin(n));
    ++n;
  }
}

void SimpleSimulator::Shed(const std::size_t step) {
  std::vector<Eigen::Vector3d> te_pos;
  std::vector<double> te_gamma;
  for (std::size_t n = 0; n < wing_pos_.num(); n++) {
    te_pos.insert(te_pos.end(), wing_pos_.iterator_at(n, wing_pos_.rows() - 1, 0),
                  wing_pos_.iterator_at(n + 1, 0, 0));
  }
  wake_pos_.prepend_row(te_pos.begin(), te_pos.end());
  if (step > 1) {
    for (std::size_t n = 0; n < wing_gamma_.num(); n++) {
      te_gamma.insert(te_gamma.end(),
                      wing_gamma_.iterator_at(n, wing_gamma_.rows() - 1, 0),
                      wing_gamma_.iterator_at(n + 1, 0, 0));
    }
    wake_gamma_.prepend_row(te_gamma.begin(), te_gamma.end());
  }
}

void SimpleSimulator::Advect(const double dt) {
  std::vector<Eigen::Vector3d> wake_vel(wake_pos_.size());

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (std::size_t i = 0; i < wake_pos_.size(); i++) {
    // wake_vel[i] = Velocity(wake_pos[i]);
    wake_vel[i] = -forward_flight_;
  }

  // for trailing edge
  std::fill(wake_vel.begin(), wake_vel.end(), -forward_flight_);
  for (std::size_t i = 0; i < wake_pos_.size(); i++) {
    wake_pos_[i] += wake_vel[i] * dt;
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
  LOG(INFO) << "Shed";
  Shed(step);
  
  const auto cpos = CollocationPoints(wing_pos_);
  const auto normal = Normals(wing_pos_);

  if (result_path_.size()) OutputPanels(step, dt);

  LOG(INFO) << "Advect";
  Advect(dt);
}

void SimpleSimulator::Run(const std::size_t steps, const double dt) {
  BuildWing();
  wake_pos_.resize(wing_pos_.num(), 0, wing_pos_.cols());
  wake_gamma_.resize(wing_gamma_.num(), 0, wing_gamma_.cols());
  for (std::size_t step = 1; step<=steps; step++) {
    MainLoop(step, dt);
  }
}

void SimpleSimulator::OutputPanels(const std::size_t step, const double dt) const {
  char filename[256];
  UVLM::proto::Snapshot2 snapshot;
  for (std::size_t n = 0; n < wing_gamma_.num(); n++) {
    UVLM::output::SimpleAppendSnapshot(
        &snapshot, wing_pos_.sheet_begin(n), wing_pos_.sheet_end(n),
        wing_gamma_.sheet_begin(n), wing_gamma_.sheet_end(n),
        wing_gamma_.cols());
    if (wake_gamma_.size()) {
      UVLM::output::SimpleAppendSnapshot(
          &snapshot, wake_pos_.sheet_begin(n), wake_pos_.sheet_end(n),
          wake_gamma_.sheet_begin(n), wake_gamma_.sheet_end(n),
          wake_gamma_.cols());
    }
  }

  // TODO change path
  sprintf(filename, "%s/%08lu", result_path_.c_str(), step);
  std::ofstream ofs(filename);
  CHECK(ofs);
  snapshot.SerializeToOstream(&ofs);
}

}  // namespace simulator
}  // namespace UVLM

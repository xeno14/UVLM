/**
 * @file simple_simulator.h
 * @brief Add description here
 */
#pragma once

#include <fstream>

#include "../../proto/uvlm.pb.h"
#include "../advect.h"
#include "../calc_load/joukowski.h"
#include "../morphing.h"
#include "../multiple_sheet/multiple_sheet.h"
#include "../recordio/recordio.h"
#include "../velocity.h"
#include "../wing/wing.h"

using multiple_sheet::MultipleSheet;

namespace UVLM {
namespace simulator {

struct WingInformation {
  std::unique_ptr<wing::WingGenerator> generator;
  Morphing morphing;
  double chord, span;
  std::size_t rows, cols;
  Eigen::Vector3d origin;
};

std::vector<Eigen::Vector3d>
    CollocationPoints(const MultipleSheet<Eigen::Vector3d>& pos);

std::vector<Eigen::Vector3d> Normals(const MultipleSheet<Eigen::Vector3d>& pos);

class SimpleSimulator {
 public:
  static const std::size_t STEP_MIN = 0;

  SimpleSimulator()
      : forward_flight_(Eigen::Vector3d::Zero()), ofs_load_(nullptr) {}
  ~SimpleSimulator() {}

  /**
   * @brief Add wing information
   *
   * @param chor
   * @todo specify wing generator
   */
  void AddWing(wing::WingGenerator* wing_generator,
               const Morphing& morphing, const double chord, const double span,
               const std::size_t rows, const std::size_t cols,
               const Eigen::Vector3d& origin);

  void Run(const std::size_t steps, const double dt);

  void set_forward_flight(const Eigen::Vector3d& v) { forward_flight_ = v; }
  void set_result_path(const std::string& path);
  void set_load_path(const std::string& loadpath);
  void set_sheet_path(const std::string& path);
  void set_advection(advect::Advection* advection) {
    advection_.reset(advection);
  }

 private:
  std::unique_ptr<advect::Advection> advection_;
  MultipleSheet<Eigen::Vector3d> wing_pos_;
  MultipleSheet<Eigen::Vector3d> wing_pos_init_;
  MultipleSheet<Eigen::Vector3d> wake_pos_;
  MultipleSheet<double> wing_gamma_;
  MultipleSheet<double> wing_gamma_prev_;
  MultipleSheet<double> wake_gamma_;
  Eigen::Vector3d forward_flight_;
  std::vector<Morphing> morphings_;
  std::vector<WingInformation> wing_info_;
  std::unique_ptr<std::ofstream> ofs_result_;
  std::unique_ptr<std::ofstream> ofs_sheet_;
  std::unique_ptr<std::ofstream> ofs_load_;
  std::unique_ptr<recordio::RecordWriter> writer_;
  std::unique_ptr<recordio::RecordWriter> sheet_writer_;

  Eigen::MatrixXd CalcMatrix(const std::vector<Eigen::Vector3d>& cpos,
                             const std::vector<Eigen::Vector3d>& normal) const;
  Eigen::VectorXd CalcRhs(const std::vector<Eigen::Vector3d>& cpos,
                          const std::vector<Eigen::Vector3d>& cpos_init,
                          const std::vector<Eigen::Vector3d>& normal,
                          const double t) const;
  void Shed(const std::size_t step);
  void Advect(const double dt);
  void CalcLoad(const std::vector<Eigen::Vector3d>& normal,
                const double t, const double dt) const;
  Eigen::Vector3d BoundVelocity(const Eigen::Vector3d& x) const {
    return UVLM::InducedVelocity(x, wing_pos_, wing_gamma_);
  }
  Eigen::Vector3d WakeVelocity(const Eigen::Vector3d& x) const {
    return UVLM::InducedVelocity(x, wake_pos_, wake_gamma_);
  }
  Eigen::Vector3d Velocity(const Eigen::Vector3d& x) const {
    return -forward_flight_ + BoundVelocity(x) + WakeVelocity(x);
  }
  void EraseOldestWake();

  void BuildWing();
  void MainLoop(const std::size_t step, const double dt);
  void OutputPanels(const std::size_t step, const double dt) const;
  void OutputSheet(const std::size_t step, const double dt) const;
  void PrepareOutputLoad();
};
}  // namespace simulator
}  // namespace UVLM

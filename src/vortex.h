
/**
 * @file vortex.h
 * @brief Add description here
 */
#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <initializer_list>
#include <vector>

using Eigen::Vector3d;

namespace UVLM {

/** @brief 点と直線の距離
 */
double DistanceLineAndPoint(const Vector3d& start,
                            const Vector3d& end,
                            const Vector3d& pos);

/** @brief Biot-Savartの法則の線分バージョン
 *
 * @param result 結果
 * @param start 始点
 * @param end 終点
 * @param pos 位置
 */
void BiotSavartLaw(Vector3d* result,
                   const Vector3d& start, const Vector3d& end,
                   const Vector3d& pos);

/** @brief Vortex filament
 * 
 * 始点・終点、循環を持つ
 */
class VortexFilament {
 public:
  VortexFilament() {}
  VortexFilament(double gamma, const Vector3d& start,
                 const Vector3d& end)
      : gamma_(gamma), start_(start), end_(end) {}
  VortexFilament(const VortexFilament& v)
      : gamma_(v.gamma_), start_(v.start_), end_(v.end_) {}

  void BiotSavartLaw(Vector3d* result, const Vector3d& pos) const;

  double gamma() const { return gamma_; }
  void set_gamma(double gamma) { gamma_ = gamma; }
  Vector3d& start() { return start_; }
  const Vector3d& start() const { return start_; }
  Vector3d& end() { return end_; }
  const Vector3d& end() const { return end_; }

 private:
  double gamma_;
  Vector3d start_, end_;
};


/** @brief 渦輪=Vortex filamentの集合
 */
class VortexRing {
 public:
  VortexRing() {}
  VortexRing(const VortexRing& v)
      : gamma_(v.gamma_), filaments_(v.filaments_), nodes_() {}

  void BiotSavartLaw(Vector3d* result, const Vector3d& pos) const;

  VortexRing& PushNode(const Vector3d& pos);
  void AssembleRing();

  double gamma() const { return gamma_; }
  void set_gamma(double gamma) { gamma_ = gamma; }
  std::vector<VortexFilament>& filaments() { return filaments_; }
  const std::vector<VortexFilament>& filaments() const { return filaments_; }

  void Clear() { filaments_.clear(); nodes_.clear(); }

 private:
  double gamma_;
  std::vector<VortexFilament> filaments_;
  std::vector<Vector3d> nodes_;
};

}  // namespace UVLM

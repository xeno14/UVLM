/**
 * @file vortex.h
 * @brief Add description here
 */
#pragma once

#include "multiple_sheet/multiple_sheet.h"
#include "vortex_kernel.h"

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <initializer_list>
#include <vector>

using Eigen::Vector3d;
using multiple_sheet::MultipleSheet;

namespace UVLM {

/** @brief 点と直線の距離
 */
double DistanceLineAndPoint(const Vector3d& start, const Vector3d& end,
                            const Vector3d& pos);

/** @brief Biot-Savartの法則の線分バージョン
 *
 * @param result 結果
 * @param start 始点
 * @param end 終点
 * @param pos 位置
 */
void BiotSavartLaw(Vector3d* result, const Vector3d& start, const Vector3d& end,
                   const Vector3d& pos);

/** @brief 渦輪=Vortex filamentの集合
 */
class VortexRing {
 public:
  static const std::size_t DEFAULT_NODE_SIZE = 4;

  VortexRing() : nodes_() {}
  VortexRing(const VortexRing& v)
      : gamma_(v.gamma_), nodes_(v.nodes_), nodes0_(v.nodes0_) {}

  /** @brief この渦がBiot-Savartの法則の法則によりつくる流れ
   *  @param[out] result 結果
   *  @param[in]  pos 流れの位置
   *  @param[in]  循環の強さ
   */
  void BiotSavartLaw(Vector3d* result, const Vector3d& pos, double gamma) const;

  void ChordwiseBiotSavartLaw(Vector3d* result, const Vector3d& pos,
                              double gamma) const;

  /** @brief この渦の循環の値を使ったBiotSavartLaw */
  inline void BiotSavartLaw(Vector3d* result, const Vector3d& pos) const {
    BiotSavartLaw(result, pos, gamma_);
  }

  inline void ChordwiseBiotSavartLaw(Vector3d* result,
                                     const Vector3d& pos) const {
    return ChordwiseBiotSavartLaw(result, pos, gamma_);
  }

  VortexRing& PushNode(const Vector3d& pos);
  VortexRing& PushNode(double x, double y, double z);
  void SaveReferenceNode();

  /** @brief 法線ベクトル
   *  @todo キャッシュする
   */
  Vector3d Normal() const;

  /** @brief tangent vector
   */
  Vector3d Tangent() const;

  /** @brief tangent vector perpendiculer to normal and tangent
   */
  Vector3d Tangent2() const;

  /** @brief 中心の位置
   *  @todo ちゃんと計算する
   */
  Vector3d Centroid() const;
  Vector3d ReferenceCentroid() const;

  double gamma() const { return gamma_; }
  void set_gamma(double gamma) { gamma_ = gamma; }
  std::vector<Vector3d>& nodes() { return nodes_; }
  const std::vector<Vector3d>& nodes() const { return nodes_; }
  const std::vector<Vector3d>& nodes0() const { return nodes0_; }

  void Clear() { nodes_.clear(); }

  /**
   * i方向(chordwise)の規格化した接ベクトルを返す
   */
  Eigen::Vector3d TanVecChord() const;

  /**
   * j方向(spanwise)の規格化した接ベクトルを返す
   */
  Eigen::Vector3d TanVecSpan() const;

  /**
   * 渦輪のchord方向の長さを求める
   */
  double CalcC() const;

  /**
   * 渦輪のspan方向の長さを求める
   */
  double CalcB() const;

  /**
   * @brief Area of ring
   * sum of triangles
   */
  double Area() const;

  /**
   * @brief Vortex impulse
   */
  Eigen::Vector3d Impulse() const;

  /**
   * @brief loop for each line segment of ring
   *
   * start and end are given for the function. the function must be formed
   * func(const Eigen::Vector3d& start, const Eigen::Vector3d& end)
   */
  inline void ForEachSegment(
      std::function<void(const Eigen::Vector3d&, const Eigen::Vector3d&)> func)
      const {
    for (std::size_t i = 0; i < nodes_.size(); i++) {
      const Vector3d& start = nodes_[(i + 1) % nodes_.size()];
      const Vector3d& end = nodes_[i];
      func(start, end);
    }
  }

  /**
   * @brief loop for each line segment of ring with reference positions
   *
   * start, end, start0 and end0 are given for the function. the function must
   * be formed
   * func(const Eigen::Vector3d& start, const Eigen::Vector3d& end)
   */
  inline void ForEachSegment(
      std::function<void(const Eigen::Vector3d&, const Eigen::Vector3d&,
                         const Eigen::Vector3d&, const Eigen::Vector3d&)> func)
      const {
    for (std::size_t i = 0; i < nodes_.size(); i++) {
      const Vector3d& start = nodes_[(i + 1) % nodes_.size()];
      const Vector3d& end = nodes_[i];
      const Vector3d& start0 = nodes0_[(i + 1) % nodes_.size()];
      const Vector3d& end0 = nodes0_[i];
      func(start, end, start0, end0);
    }
  }

  /**
   * angle of attack
   */
  double AngleOfAttack(const Eigen::Vector3d& Q) const;

 private:
  double gamma_;
  std::vector<Vector3d> nodes_;
  std::vector<Vector3d> nodes0_;  // reference nodes (immutable)
};

inline Eigen::Vector3d VORTEX(const Eigen::Vector3d& x,
                              const Eigen::Vector3d& x1,
                              const Eigen::Vector3d& x2, double gamma) {
  // (10.16)
  // impl p. 584
  Eigen::Vector3d res;
  UVLM::BiotSavartLaw(&res, x1, x2, x);
  res *= gamma;
  return res;
}

Eigen::Vector3d VORING(const Eigen::Vector3d& x,
                       const MultipleSheet<Eigen::Vector3d>& pos,
                       const MultipleSheet<double>& gamma, std::size_t n,
                       std::size_t i, std::size_t j);

Eigen::Vector3d VORING(const vortex_kernel::VortexKernel& kernel,
                       const Eigen::Vector3d& x,
                       const MultipleSheet<Eigen::Vector3d>& pos,
                       const MultipleSheet<double>& gamma, std::size_t n,
                       std::size_t i, std::size_t j);

}  // namespace UVLM

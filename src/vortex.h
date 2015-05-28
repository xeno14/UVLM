
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


/** @brief 渦輪=Vortex filamentの集合
 */
class VortexRing {
 public:
  VortexRing() {}
  VortexRing(const VortexRing& v)
      : gamma_(v.gamma_), nodes_(v.nodes_), nodes0_(v.nodes0_) {}

  /** @brief この渦がBiot-Savartの法則の法則によりつくる流れ
   *  @param[out] result 結果
   *  @param[in]  pos 流れの位置
   *  @param[in]  循環の強さ
   */
  void BiotSavartLaw(Vector3d* result, const Vector3d& pos, double gamma) const;

  /** @brief この渦の循環の値を使ったBiotSavartLaw */
  inline void BiotSavartLaw(Vector3d* result, const Vector3d& pos) const {
    BiotSavartLaw(result, pos, gamma_);
  }

  VortexRing& PushNode(const Vector3d& pos);
  VortexRing& PushNode(double x, double y, double z);
  void SaveReferenceNode();

  /** @brief 法線ベクトル */
  Vector3d Normal() const;
  /** @brief 中心の位置 
   *  @todo ちゃんと計算する
   */
  Vector3d Centroid() const;

  double gamma() const { return gamma_; }
  void set_gamma(double gamma) { gamma_ = gamma; }
  std::vector<Vector3d>& nodes() { return nodes_; }
  const std::vector<Vector3d>& nodes() const { return nodes_; }
  const std::vector<Vector3d>& nodes0() const { return nodes0_; }

  void Clear() { nodes_.clear(); }

 private:
  double gamma_;
  std::vector<Vector3d> nodes_;
  std::vector<Vector3d> nodes0_; // reference nodes (immutable)
};

}  // namespace UVLM

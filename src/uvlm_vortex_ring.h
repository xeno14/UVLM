/**
 * @file uvlm_vortex_ring.h
 * @brief Add description here
 */
#pragma once

#include "vortex.h"

#include <vector>

namespace UVLM {
namespace internal {

/** @brief 渦輪によってつくられる流れ
 *  @param[out] result 流速
 *  @param[in]  pos 位置
 *  @param[in]  vortices 渦輪
 */
void InducedVelocityByVortices(Eigen::Vector3d* const result,
                               const Eigen::Vector3d& pos,
                               const std::vector<VortexRing>& vortices);

}  // internal

/** @brief
 *
 *  こんなことができます
 *  - ある位置における流速の計算
 *  - 渦の世代管理
 *  - Trailing edgeに関する処理
 *  - Wing tipに関する処理
 */
class UVLMVortexRing {
 public:
  UVLMVortexRing() {}
  ~UVLMVortexRing() = default;

  /** @brief 翼上の点からbound vortexの節の位置を決める
   *  @pre 右半分(y>0)で、左上(x,y=0,0)から右下(x,y=S,C)に走査する並びで、yが動く
   */
  void InitWing(const std::vector<Eigen::Vector3d>& pos, std::size_t cols); 

  /** @brief 全ての渦によりつくられる流れ */
  void InducedVelocity(Eigen::Vector3d* const result,
                       const Eigen::Vector3d& pos) const;
  /** @brief 翼上の渦によりつくられる流れ */
  inline void InducedVelocityByBound(Eigen::Vector3d* const result,
                                     const Eigen::Vector3d& pos) const {
    internal::InducedVelocityByVortices(result, pos, bound_vortices_);
  }
  /** @brief 後流の渦によりつくられる流れ */
  inline void InducedVelocityByWake(Eigen::Vector3d* const result,
                                    const Eigen::Vector3d& pos) const {
    internal::InducedVelocityByVortices(result, pos, wake_vortices_);
  }

  std::vector<VortexRing>& bound_vortices() { return bound_vortices_; }
  std::vector<VortexRing>& wake_vortices() { return wake_vortices_; }
  const std::vector<VortexRing>& bound_vortices() const {
    return bound_vortices_;
  }
  const std::vector<VortexRing>& wake_vortices() const { return wake_vortices_; }

  std::size_t rows() const { return rows_; }
  std::size_t cols() const { return cols_; }

 private:
  std::vector<VortexRing> bound_vortices_;
  std::vector<VortexRing> wake_vortices_;
  std::size_t cols_, rows_;
};

}  // namespace UVLM

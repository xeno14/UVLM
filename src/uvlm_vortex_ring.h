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
  UVLMVortexRing(std::size_t cols, std::size_t rows) : cols_(cols), rows_(rows) {}
  ~UVLMVortexRing() = default;

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
  const std::size_t cols_, rows_;
};

}  // namespace UVLM

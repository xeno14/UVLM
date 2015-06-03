/**
 * @file uvlm_vortex_ring.h
 * @brief Add description here
 */
#pragma once

#include "vortex.h"
#include "iterator.h"

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

inline void CopySymmetry(Eigen::Vector3d* const to,
                         const Eigen::Vector3d& from) {
  *to << from.x(), -from.y(), from.z();
}

}  // internal

/** @brief
 *
 *  @todo 位置は原点からの相対的な位置にする？
 *
 *  こんなことができます
 *  - ある位置における流速の計算
 *  - 渦の世代管理
 *  - Trailing edgeに関する処理
 *  - Wing tipに関する処理
 */
class UVLMVortexRing {
 public:
  typedef Iterator<std::vector<VortexRing>> iterator;

  UVLMVortexRing()
      : bound_vortices_(),
        wake_vortices_(),
        cols_(0),
        rows_(0),
        origin_(Eigen::Vector3d::Zero()) {}
  ~UVLMVortexRing() = default;

  /** @brief 翼上の点からbound vortexの節の位置を決める
   *
   *  注意：cols は渦の個数を表す。なので、posは(cols+1)*(rows+1)個の
   *  要素数でないといけない。
   *  @param pos 位置
   *  @param cols 渦の個数
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

  /** @brief bound vorticesの位置を原点を通るx-z平面に対し面対称にする */
  void PlaneSymmetry();

  std::pair<iterator, iterator> TrailingEdgeIterators() {
    std::size_t first = 0 + (rows_ - 1) * 2 * cols_;
    std::size_t last =  0 + rows_ * 2 * cols_;
    return std::make_pair(iterator(bound_vortices_, first, 1),
                          iterator(bound_vortices_, last, 1));
  }

  std::vector<VortexRing>& bound_vortices() { return bound_vortices_; }
  const std::vector<VortexRing>& bound_vortices() const { return bound_vortices_; }

  std::vector<VortexRing>& wake_vortices() { return wake_vortices_; }
  const std::vector<VortexRing>& wake_vortices() const { return wake_vortices_; }

  std::size_t rows() const { return rows_; }
  std::size_t cols() const { return cols_; }

 private:
  std::vector<VortexRing> bound_vortices_;
  std::vector<VortexRing> wake_vortices_;
  std::size_t cols_, rows_;
  Eigen::Vector3d origin_;
};

}  // namespace UVLM

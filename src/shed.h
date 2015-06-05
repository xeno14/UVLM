
/**
 * @file shed.h
 * @brief Add description here
 */
#pragma once

#include "vortex.h"
#include "morphing.h"
#include "uvlm_vortex_ring.h"

namespace UVLM {
namespace internal {

/** @brief Euler-methodで移流させる
 *  @param[out] target 位置
 *  @param[in]  vel 速度
 *  @param[in]  dt  時間刻み
 */
void Advect(Eigen::Vector3d* target, const Eigen::Vector3d& vel,
            const double dt);

/** @biref Trailing edgeの渦を一つ放出する
 *  @param[out] result 放出された渦
 *  @param[in]  target 放出される前のtrailing edgeにある渦
 *  @param[in]  rings 渦
 *  @param[in]  dt 時間刻み
 */
void ShedSingleAtTrailingEdge(VortexRing* result, const VortexRing& target,
                              const UVLMVortexRing& rings,
                              const Morphing& morphing,
                              const Eigen::Vector3d& Vinfty, const double t,
                              const double dt);

/** @brief 移流の実装
 *
 *  result に移流した結果を入れる
 *  @todo ここが律速になるはずなので並列化する
 *  @param[out] result 移流の結果
 *  @param[in]  rings 渦
 *  @param[in]  dt 時間刻み
 */
void AdvectWakeImpl(std::vector<VortexRing>* result,
                    const UVLMVortexRing& rings, const Eigen::Vector3d& Vinfty,
                    const double dt);

}  // namespace internal

template <class InputIterator, class OutputIterator>
void ShedAtTrailingEdge(InputIterator first, InputIterator last,
                        OutputIterator result, const UVLMVortexRing& rings,
                        const Morphing& morphing,
                        const Eigen::Vector3d& Vinfty, const double t, const double dt) {
  while (first != last) {
    internal::ShedSingleAtTrailingEdge(&(*result), *first, rings, morphing, Vinfty, t, dt);
    ++first; ++result;
  }
}

void AdvectWake(UVLMVortexRing* rings, const Eigen::Vector3d& Vinfty,
                const double dt);

}  // namespace UVLM


/**
 * @file shed.h
 * @brief Add description here
 */
#pragma once

#include "vortex.h"
#include "uvlm_vortex_ring.h"

namespace UVLM {
namespace internal {

/** @brief Euler-methodで移流させる
 *  @param[out] target 位置
 *  @param[in]  vel 速度
 *  @param[in]  dt  時間刻み
 */
inline void Advect(Eigen::Vector3d* const target, const Eigen::Vector3d& vel,
            const double dt) {
  *target = *target + vel * dt;
}

/** @biref Trailing edgeの渦を一つ放出する
 *  @param[out] result 放出された渦
 *  @param[in]  target 放出される前のtrailing edgeにある渦
 *  @param[in]  rings 渦
 *  @param[in]  dt 時間刻み
 */
void ShedSingleAtTrailingEdge(VortexRing* result, const VortexRing& target,
                              const UVLMVortexRing& rings,
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

template <class InputIterator, class OutputIterator>
void AdvectWakeImpl(OutputIterator wake_first, OutputIterator wake_last,
                    InputIterator vortices_first, InputIterator vortices_last,
                    const Eigen::Vector3d& Vinfty, const double dt);

}  // namespace internal

template <class InputIterator, class OutputIterator>
void ShedAtTrailingEdge(InputIterator first, InputIterator last,
                        OutputIterator result, const UVLMVortexRing& rings,
                        const Eigen::Vector3d& Vinfty, const double t, const double dt) {
  while (first != last) {
    internal::ShedSingleAtTrailingEdge(&(*result), *first, rings, Vinfty, t, dt);
    ++first; ++result;
  }
}

template <class InputIterator1, class InputIterator2, class OutputIterator>
void ShedAtTrailingEdge(InputIterator1 edge_first, InputIterator1 edge_last,
                        OutputIterator result, InputIterator2 vortices_first,
                        InputIterator2 vortices_last,
                        const Eigen::Vector3d& Vinfty, const double t,
                        const double dt);

template <class InputIterator, class OutputIterator>
void AttachShedVorticesToEdge(InputIterator edge_first, InputIterator edge_last,
                              OutputIterator result);

void AdvectWake(UVLMVortexRing* rings, const Eigen::Vector3d& Vinfty,
                const double dt);

template <class InputIterator, class OutputIterator>
void AdvectWake(OutputIterator wake_first, OutputIterator wake_last,
                InputIterator vortices_first, InputIterator vortices_last,
                const Eigen::Vector3d& Vinfty, const double dt);

template <class InputIterator>
void InducedVelocity(Eigen::Vector3d* const result,
                     const Eigen::Vector3d& pos,
                     InputIterator first, InputIterator last);

}  // namespace UVLM

#include "shed.inl"

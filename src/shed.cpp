/**
 * @file shed.cpp
 * @brief Add description here
 */

#include "shed.h"


namespace UVLM {
namespace internal {

/** @brief Euler-methodで移流させる
 *  @param[out] target 位置
 *  @param[in]  vel 速度
 *  @param[in]  dt  時間刻み
 */
void Advect(Eigen::Vector3d* target, const Eigen::Vector3d& vel,
            const double dt) {
  *target += vel* dt;
}

/** @biref Trailing edgeの渦を一つ放出する
 *  @param[out] result 放出された渦
 *  @param[in]  target 放出される前のtrailing edgeにある渦
 *  @param[in]  rings 渦
 *  @param[in]  dt 時間刻み
 */
void ShedSingleAtTrailingEdge(VortexRing* result, const VortexRing& target,
                              const UVLMVortexRing& rings,
                              const Eigen::Vector3d& Vinfty, const double dt) {
  // before
  // 3--2=3'---2'
  // |   |     |
  // 0--1=0'---1'
  //     after
  result->nodes().resize(VortexRing::DEFAULT_NODE_SIZE);
  result->nodes()[0] = target.nodes()[1];
  result->nodes()[3] = target.nodes()[2];
  result->nodes()[1] = target.nodes()[1];  // 後で移流される
  result->nodes()[2] = target.nodes()[2];  // 後で移流される

  Eigen::Vector3d v1, v2;
  rings.InducedVelocity(&v1, result->nodes()[1]);
  rings.InducedVelocity(&v2, result->nodes()[2]);
  v1 += Vinfty;
  v2 += Vinfty;
  Advect(&result->nodes()[1], v1, dt);
  Advect(&result->nodes()[2], v2, dt);
}

}  // namespace internal

void AdvectWake(UVLMVortexRing* rings, const Eigen::Vector3d& Vinfty,
                const double dt) {
  for (auto& wake : rings->wake_vortices()) {
    Eigen::Vector3d velocity;
    for (auto& node : wake.nodes()) {
      rings->InducedVelocity(&velocity, node);
      velocity += Vinfty;
      internal::Advect(&node, velocity, dt);
    }
  }
}

}  // namespace UVLM

#pragma once

#include <iostream>

namespace UVLM {
namespace internal {

template <class InputIterator, class OutputIterator>
void AdvectWakeImpl(OutputIterator wake_first, OutputIterator wake_last,
                    InputIterator vortices_first, InputIterator vortices_last,
                    const Eigen::Vector3d& Vinfty, const double dt) {
  for (auto wake = wake_first; wake != wake_last; ++wake) {
    Eigen::Vector3d velocity;
    for (auto& node : wake->nodes()) {
      InducedVelocity(&velocity, node, vortices_first, vortices_last);
      velocity += Vinfty;
      internal::Advect(&node, velocity, dt);
    }
  }
}

}  // namespace internal

template <class InputIterator, class OutputIterator>
void AdvectWake(OutputIterator wake_first, OutputIterator wake_last,
                InputIterator vortices_first, InputIterator vortices_last,
                const UVLMVortexRing& rings,
                const Eigen::Vector3d& Vinfty, const double dt) {
  std::vector<VortexRing> new_wake(wake_first, wake_last);
  // TODO こいつが悪い
  // internal::AdvectWakeImpl(new_wake.begin(), new_wake.end(), vortices_first,
  //                          vortices_last, Vinfty, dt);
  internal::AdvectWakeImpl(&new_wake, rings, Vinfty, dt);
  std::copy(new_wake.begin(), new_wake.end(), wake_first);
}

template <class InputIterator>
void InducedVelocity(Eigen::Vector3d* const result,
                     const Eigen::Vector3d& pos,
                     InputIterator first, InputIterator last) {
  *result = Eigen::Vector3d::Zero();
  Eigen::Vector3d tmp;
  while (first != last) {
    first->BiotSavartLaw(&tmp, pos);
    *result += tmp;
    ++first;
  }
}

template <class InputIterator1, class InputIterator2, class OutputIterator>
void ShedAtTrailingEdge(InputIterator1 edge_first, InputIterator1 edge_last,
                        OutputIterator result, InputIterator2 vortices_first,
                        InputIterator2 vortices_last,
                        const UVLMVortexRing& rings,
                        const Eigen::Vector3d& Vinfty, const double t,
                        const double dt) {
  for (auto edge = edge_first; edge != edge_last; ++edge, ++result) {
    *result = *edge;
    Eigen::Vector3d velocity;
    // 移流する
    for (auto& node : result->nodes()) {
      // TODO ここがおかしい
      // InducedVelocity(&velocity, node, vortices_first, vortices_last);
      rings.InducedVelocity(&velocity, node);

      velocity += Vinfty;
      internal::Advect(&node, velocity, dt);
    }
  }
}

template <class InputIterator, class OutputIterator>
void AttachShedVorticesToEdge(InputIterator edge_first, InputIterator edge_last,
                              OutputIterator result) {
  while (edge_first != edge_last) {
    result->nodes()[0] = edge_first->nodes()[1];
    result->nodes()[3] = edge_first->nodes()[2];
    ++edge_first; ++result;
  }
}

}  // namespace UVLM

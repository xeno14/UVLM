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

template <class InputIterator>
void AdvectKernel(VortexRing* v, InputIterator first, InputIterator last,
    const Eigen::Vector3d& freestream, const double dt) {
  Eigen::Vector3d vel;
  for (auto& node : v->nodes()) {
    Eigen::Vector3d vel;
    InducedVelocity(&vel, node, first, last);
    vel += freestream;
    node += vel * dt;
  }
}

}  // namespace internal

template <class InputIterator>
void AdvectWake(std::vector<UVLM::VortexRing>* wake,
                InputIterator vortices_first, InputIterator vortices_last,
                const UVLMVortexRing& rings,
                const Eigen::Vector3d& Vinfty, const double dt) {
  // TODO こいつが悪い
  // internal::AdvectWakeImpl(new_wake.begin(), new_wake.end(), vortices_first,
  //                          vortices_last, Vinfty, dt);
  internal::AdvectWakeImpl(wake, rings, Vinfty, dt);
}

template <class InputIterator>
void InducedVelocity(Eigen::Vector3d* result,
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

template <class InputIterator, class OutputIterator>
void Advect(InputIterator vortex_first, InputIterator vortex_last,
            OutputIterator wake_first, OutputIterator wake_last,
            const Eigen::Vector3d& freestream, const double dt) {
  for (auto wake=wake_first; wake != wake_last; ++wake) {
    internal::AdvectKernel(&(*wake), vortex_first, vortex_last, freestream, dt);
  }
}

template <class InputIterator, class OutputIterator>
void AdvectParallel(InputIterator vortex_first, InputIterator vortex_last,
            OutputIterator wake_first, OutputIterator wake_last,
            const Eigen::Vector3d& freestream, const double dt) {
  const std::size_t num = std::distance(wake_first, wake_last);
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (std::size_t i = 0; i < num; i++) {
    auto wake = wake_first + i;
    internal::AdvectKernel(&(*wake), vortex_first, vortex_last, freestream, dt);
  }
}

template <class InputIterator>
void ChordwiseInducedVelocity(Eigen::Vector3d* const result,
                              const Eigen::Vector3d& pos, InputIterator first,
                              InputIterator last) {
  *result = Eigen::Vector3d::Zero();
  Eigen::Vector3d tmp;
  while (first != last) {
    first->ChordwiseBiotSavartLaw(&tmp, pos);
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
      // rings.InducedVelocity(&velocity, node);

      // velocity = Vinfty;
      velocity = Eigen::Vector3d(1, 0, 0);
      internal::Advect(&node, velocity, dt);
    }
  }
}

template <class InputIterator, class OutputIterator>
void ConnectTrailingEdge(InputIterator edge_first, InputIterator edge_last,
                         OutputIterator wake) {
  for (auto edge = edge_first; edge != edge_last; ++edge, ++wake) {
    wake->nodes()[0] = edge->nodes()[1];
    wake->nodes()[3] = edge->nodes()[2];
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

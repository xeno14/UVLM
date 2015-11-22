#pragma once

#include <iostream>

namespace UVLM {
namespace internal {

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

template <class InputIterator, class OutputIterator>
void ConnectTrailingEdge(InputIterator edge_first, InputIterator edge_last,
                         OutputIterator wake) {
  for (auto edge = edge_first; edge != edge_last; ++edge, ++wake) {
    wake->nodes()[0] = edge->nodes()[1];
    wake->nodes()[3] = edge->nodes()[2];
  }
}

}  // namespace UVLM

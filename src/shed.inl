#pragma once

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
                const Eigen::Vector3d& Vinfty, const double dt) {
  std::vector<VortexRing> new_wake(wake_first, wake_last);
  internal::AdvectWakeImpl(new_wake.begin(), new_wake.end(), vortices_first,
                           vortices_last, Vinfty, dt);
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

}  // namespace UVLM

#pragma once

namespace UVLM {
namespace internal {

template <class InputIterator, class OutputIterator>
void AdvectWakeImpl(OutputIterator wake_first, OutputIterator wake_last,
                    InputIterator vortices_first, InputIterator vortices_last,
                    const Eigen::Vector3d& Vinfty, const double dt) {
  while (wake_first != wake_last) {
    Eigen::Vector3d velocity;
    for (auto& node : wake_first->nodes()) {
      InducedVelocity(&velocity, node, vortices_first, vortices_last);
      velocity += Vinfty;
      internal::Advect(&node, velocity, dt);
    }
    ++wake_first;
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

}  // namespace UVLM

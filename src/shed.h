
/**
 * @file shed.h
 * @brief Add description here
 */
#pragma once
#include "vortex.h"
#include "uvlm_vortex_ring.h"
#include "vortex_container.h"

namespace UVLM {
namespace internal {

template <class InputIterator>
void AdvectKernel(VortexRing* v, InputIterator first, InputIterator last,
    const Eigen::Vector3d& freestream, const double dt);

}  // namespace internal

template <class OutputIterator>
void ConnectTrailingEdge(OutputIterator wake, const VortexContainer& c);

template <class InputIterator, class OutputIterator>
void Advect(InputIterator vortex_first, InputIterator vortex_last,
            OutputIterator wake_first, OutputIterator wake_last,
            const Eigen::Vector3d& freestream, const double dt);

template <class InputIterator, class OutputIterator>
void AdvectParallel(InputIterator vortex_first, InputIterator vortex_last,
                    OutputIterator wake_first, OutputIterator wake_last,
                    const Eigen::Vector3d& freestream, const double dt);

template <class InputIterator>
void InducedVelocity(Eigen::Vector3d* const result,
                     const Eigen::Vector3d& pos,
                     InputIterator first, InputIterator last);

template <class InputIterator>
void ChordwiseInducedVelocity(Eigen::Vector3d* const result,
                              const Eigen::Vector3d& pos, InputIterator first,
                              InputIterator last);

}  // namespace UVLM

#include "shed.inl"

/**
 * @file output.cpp
 * @brief Add description here
 */

#include "output.h"

namespace UVLM {
namespace output {

void SheetToSnapshot(UVLM::proto::Snapshot2* snapshot,
                     const MultipleSheet<Eigen::Vector3d>& pos,
                     const MultipleSheet<double>& gamma) {
  for (const auto& index : gamma.list_index()) {
    std::size_t n, i, j;
    std::tie(std::ignore, n, i, j) = index;

    // TODO wrap here
    const auto& p0 = pos.at(n, i, j);
    const auto& p1 = pos.at(n, i, j + 1);
    const auto& p2 = pos.at(n, i + 1, j + 1);
    const auto& p3 = pos.at(n, i + 1, j);
    const double g = gamma.at(n, i, j);

    UVLM::VortexRing v;
    // TODO why opposite?
    v.PushNode(p3).PushNode(p2).PushNode(p1).PushNode(p0);
    // v.PushNode(p0).PushNode(p1).PushNode(p2).PushNode(p3);
    v.set_gamma(g);
    snapshot->add_vortices()->CopyFrom(::UVLM::VortexRingToProto(v));
  }
}

}  // namespace output
}  // namespace UVLM

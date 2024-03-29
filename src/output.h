
/**
 * @file output.h
 * @brief Add description here
 */
#pragma once

#include "../proto/uvlm.pb.h"
#include "multiple_sheet/multiple_sheet.h"
#include "proto_adaptor.h"

#include <glog/logging.h>
#include <fstream>
#include <iterator>
#include <type_traits>

using multiple_sheet::MultipleSheet;

namespace UVLM {
namespace output {

/**
 * @brief Append list elements to snapshot
 * @todo move to adaptor?
 */
template <class InputIteraotr1, class InputIteraotr2>
void SimpleAppendSnapshot(UVLM::proto::Snapshot2* snapshot,
                      InputIteraotr1 pos_first, InputIteraotr1 pos_last,
                      InputIteraotr2 gamma_first, InputIteraotr2 gamma_last,
                      const std::size_t cols) {
  // TODO type check
  // ITERATOR_VALUETYPE_CHECK(Eigen::Vector3d, InputIteraotr1);
  // ITERATOR_VALUETYPE_CHECK(double, InputIteraotr2);

  const std::size_t pos_size = std::distance(pos_first, pos_last);
  const std::size_t gamma_size = std::distance(gamma_first, gamma_last);
  CHECK(pos_size == (cols + 1) * (gamma_size / cols + 1)) << " "
    "Invalid size: pos size=" << pos_size << " "
    "gamma size=" << gamma_size;

  for (std::size_t K = 0; K < gamma_size; ++K) {
    auto gamma = gamma_first + K;
    std::size_t i = K / cols;
    std::size_t j = K % cols;
    auto pos = pos_first + j + i * (cols + 1);
    UVLM::VortexRing v;
    v.PushNode(*(pos))
        .PushNode(*(pos + 1))
        .PushNode(*(pos + 1 + cols + 1))
        .PushNode(*(pos + cols + 1));
    v.set_gamma(*gamma);
    snapshot->add_vortices()->CopyFrom(::UVLM::VortexRingToProto(v));
  }
}

void SheetToSnapshot(UVLM::proto::Snapshot2* snapshot,
                            const MultipleSheet<Eigen::Vector3d>& pos,
                            const MultipleSheet<double>& gamma);

}  // namespace output
}  // namespace UVLM

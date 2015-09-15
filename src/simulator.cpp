/**
 * @file simulator.cpp
 * @brief Add description here
 */

#include "simulator.h"

#include "../proto/uvlm.pb.h"
#include "uvlm_vortex_ring.h"
#include "morphing.h"
#include "linear.h"
#include "shed.h"
#include "proto_adaptor.h"
#include "util.h"
#include "vortex_container.h"
#include "wing_builder.h"

#include <fstream>
#include <glog/logging.h>
#include <vector>
#include <utility>


namespace {

std::vector<UVLM::proto::Wing> wings;
std::vector<UVLM::Morphing> morphings;
auto vortices = std::make_shared<std::vector<UVLM::VortexRing>>();
std::vector<UVLM::VortexContainer> containers;
UVLM::UVLMVortexRing rings;
Eigen::Vector3d inlet(1, 0, 0);
std::string output_path;

}  // anonymous namespace

namespace UVLM {
namespace simulator {
namespace internal {

std::size_t WakeOffset() {
  static std::size_t wake_offset = 0;
  if (wake_offset) return wake_offset;
  return wake_offset =
             ::UVLM::CountTotalSize(containers.begin(), containers.end());
}

void CheckReady() {
  if (morphings.size() == 0) LOG(ERROR) << "No morphing";
  if (wings.size() == 0 || containers.size() == 0) LOG(ERROR) << "No wing";
}

void CreateContainers() {
  UVLM::WingBuilder builder(&containers, vortices);
  for (const auto& wing : wings) {
    builder.AddWing(wing);
  }
  builder.Build();
} 

auto ShedProcess(const double dt) {
  // shed[i]: 翼iから放出される渦
  std::vector<std::vector<UVLM::VortexRing>> shed(containers.size());

  // 放出する渦の計算
  for (std::size_t i = 0; i < containers.size(); i++) {
    auto edge_first = containers[i].edge_begin();
    auto edge_last = containers[i].edge_end();
    shed[i].resize(std::distance(edge_first, edge_last));
        UVLM::ShedAtTrailingEdge(edge_first, edge_last, shed[i].begin(),
                                 vortices->cbegin(), vortices->cend(), rings,
                                 inlet, 0, dt);
  }
  return shed;
}

void MorphingProcess(const double t, const double dt) {
  for (std::size_t i = 0; i < containers.size(); i++) {
    for (auto& vortex : containers[i]) {
      for (std::size_t j = 0; j < vortex.nodes().size(); j++) {
        morphings[i].Perfome(&vortex.nodes()[j], vortex.nodes0()[j], t, dt);
      }
    }
  }
}

void AdvectProcess(const double dt) {
  const std::size_t wake_offset = WakeOffset(); 
  UVLM::AdvectWake(vortices->begin() + wake_offset, vortices->end(),
      vortices->cbegin(), vortices->cend(), rings, inlet, dt);
}

void AppendShedProcess(std::vector<std::vector<UVLM::VortexRing>>* shed) {
  for (std::size_t i = 0; i < containers.size(); i++) {
    UVLM::AttachShedVorticesToEdge(containers[i].edge_begin(),
        containers[i].edge_end(),
        shed->at(i).begin());
  }
  for (std::size_t i = 0; i < shed->size(); i++) {
    vortices->insert(vortices->end(), shed->at(i).cbegin(), shed->at(i).cend());
  }
}

void SolveLinearProblem(double t) {
  const std::size_t wake_offset = internal::WakeOffset();
  auto morphing = *(morphings.begin());
  auto gamma = ::UVLM::SolveLinearProblem(
      vortices->begin(), vortices->begin() + wake_offset,
      vortices->cbegin() + wake_offset, vortices->cend(), inlet, morphing, t);
  for (std::size_t i = 0; i < wake_offset; i++) {
    vortices->at(i).set_gamma(gamma(i));
  }
}

void OutputSnapshot2(const std::size_t step, const double t) {
  UVLM::proto::Snapshot2 snapshot;
  snapshot.set_t(t);

  for (const auto& container : containers) {
    auto* shape = snapshot.add_container_shapes();
    shape->set_rows(container.rows());
    shape->set_cols(container.cols());
    shape->set_id(container.id());
    // shape->set_origin(Vector3dToPoint(container->origin()));
  }
  for (const auto& vortex : *vortices) {
    snapshot.add_vortices()->CopyFrom(::UVLM::VortexRingToProto(vortex));
  }
  char filename[256];
  sprintf(filename, "%s/%08lu", output_path.c_str(), step); 
  std::ofstream ofs(filename, std::ios::binary);
  CHECK((bool)ofs) << "Unable to open " << filename;
  snapshot.SerializeToOstream(&ofs);
}

}  // namespace internal

void InitSimulator() {
  // do nothing
}

void AddWing(const ::UVLM::proto::Wing& wing, const UVLM::Morphing& m) {
  wings.emplace_back(wing);
  morphings.emplace_back(m);
}

void SetInlet(double x, double y, double z) {
  inlet = Eigen::Vector3d(x, y, z);
}

void SetOutputPath(const std::string& path) {
  output_path = path;
}

void Start(const std::size_t steps, const double dt) {
  internal::CheckReady();
  internal::CreateContainers();
  rings.bound_vortices() = *vortices;

  // Main loop
  for (std::size_t step=0; step<steps; ++step) {
    const double t = dt * step;
    const std::size_t wake_offset = internal::WakeOffset();
    // 最初のステップは放出と移流を行わない
    // Wake vortices process
    if (vortices->begin() + wake_offset == vortices->end()) {
      LOG(INFO) << "Shed";
      auto shed = internal::ShedProcess(dt);

      LOG(INFO) << "Advect";
      internal::AdvectProcess(dt);

      LOG(INFO) << "Morphing";
      internal::MorphingProcess(t, dt);

      LOG(INFO) << "Append Shed";
      internal::AppendShedProcess(&shed);
    }
    // Bound vortices process
    LOG(INFO) << "Linear problem";
    internal::SolveLinearProblem(t);

    LOG(INFO) << "Output";
    internal::OutputSnapshot2(step, t);

    // TODO remove
    std::copy(vortices->begin(), vortices->begin() + wake_offset,
              rings.bound_vortices().begin());
  }
}

}  // namespace simulator
}  // namespace UVLM

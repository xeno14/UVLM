/**
 * @file simulator.cpp
 * @brief Add description here
 */

#include "simulator.h"

#include "../proto/uvlm.pb.h"
#include "calc_load/calc_load.h"
#include "uvlm_vortex_ring.h"
#include "morphing.h"
#include "linear.h"
#include "shed.h"
#include "proto_adaptor.h"
#include "util.h"
#include "vortex_container.h"
#include "wing_builder.h"

#include <fstream>
#include <gflags/gflags.h>
#include <glog/logging.h>
#include <sstream>
#include <vector>
#include <utility>

DEFINE_string(run_name, "", "name of run");
DEFINE_bool(disable_wake, false, "disable wake shedding");

namespace {

std::vector<UVLM::proto::Wing> wings;
std::vector<UVLM::Morphing> morphings;
auto vortices = std::make_shared<std::vector<UVLM::VortexRing>>();
std::vector<UVLM::VortexContainer> containers;
std::vector<UVLM::VortexContainer> containers_prev;
UVLM::UVLMVortexRing rings;
Eigen::Vector3d inlet(1, 0, 0);
std::string output_path;
std::string output_load_path;

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
  if (morphings.size() == 0) LOG(FATAL) << FLAGS_run_name << " " << "No morphing";
  if (wings.size() == 0) LOG(FATAL) << FLAGS_run_name << " " << "No wing";

  std::ofstream ofs_load(output_load_path);
  // TODO check output path is writable
}

void CreateContainers() {
  UVLM::WingBuilder builder(&containers, vortices);
  for (const auto& wing : wings) {
    builder.AddWing(wing);
  }
  builder.Build();

  for (std::size_t i = 0; i < containers.size(); i++)
    LOG(INFO) << FLAGS_run_name << " container[" << i
              << "]: " << containers[i].rows() << " x " << containers[i].cols();
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
                             vortices->cbegin(), vortices->cend(), rings, inlet,
                             0, dt);
  }
  return shed;
}

void MorphingProcess(const double t) {
  for (std::size_t i = 0; i < containers.size(); i++) {
    for (auto& vortex : containers[i]) {
      for (std::size_t j = 0; j < vortex.nodes().size(); j++) {
        morphings[i].Perfome(&vortex.nodes()[j], vortex.nodes0()[j], t);
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

/** 
 * @todo multiple morphings
 */
void SolveLinearProblem(double t) {
  const std::size_t wake_offset = internal::WakeOffset();
  auto morphing = *(morphings.begin());
  auto gamma = ::UVLM::SolveLinearProblem(
      vortices->cbegin(), vortices->cbegin() + wake_offset,
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
    shape->set_chord(container.chord());
    shape->set_span(container.span());
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

void CalcLoadProcess(const double t, const double dt) {
  std::vector<Eigen::Vector3d> loads;
  std::stringstream line;
  for (std::size_t i=0; i<containers.size(); i++) {
    const auto& c = containers[i];
    const auto& c_prev = containers_prev[i];
    const auto& m = morphings[i];
    const double rho = 1;
    const double S = c.chord() * c.span();
    auto wake = GetWake(containers);
    auto load = UVLM::calc_load::CalcLoad(c, c_prev, wake.first, wake.second, m,
                                          inlet, rho, t, dt);
    const double U = inlet.norm();
    auto coeff = load.F / (0.5 * rho * U * U * S);

    line << t << "\t" << coeff.x() << "\t" << coeff.y() << "\t" << coeff.z()
      << "\t" << load.Pin << "\t" << load.Pout;
  }
  std::ofstream ofs(output_load_path, std::ios::app);
  CHECK((bool)ofs) << "Unable to open " << output_load_path;
  ofs << line.str() << std::endl;
}

}  // namespace internal

void InitSimulator() {
  // do nothing
}

void AddWing(const ::UVLM::proto::Wing& wing, const ::UVLM::Morphing& m) {
  wings.emplace_back(wing);
  morphings.emplace_back(m);
}

void SetInlet(double x, double y, double z) {
  inlet = Eigen::Vector3d(x, y, z);
}

void SetOutputPath(const std::string& path) { output_path = path; }

void SetOutputLoadPath(const std::string& path) { output_load_path = path; }

void Start(const std::size_t steps, const double dt) {
  internal::CheckReady();
  internal::CreateContainers();
  rings.bound_vortices() = *vortices;

  internal::MorphingProcess(0);

  // Main loop
  for (std::size_t step = 1; step <= steps; ++step) {
    LOG(INFO) << FLAGS_run_name << " " << "step: " << step;
    const double t = dt * step;
    const std::size_t wake_offset = internal::WakeOffset();

    // Save circulations of bound vortices at the previous step
    containers_prev.resize(containers.size());
    CopyContainers(containers.begin(), containers.end(),
                   containers_prev.begin());

    // 連立方程式を解いて翼の上の循環を求める
    internal::SolveLinearProblem(t);

    // TODO remove rings
    std::copy(vortices->begin(), vortices->begin() + wake_offset,
              rings.bound_vortices().begin());
    internal::OutputSnapshot2(step, t);
    if (output_load_path.size()) {
      internal::CalcLoadProcess(t, dt);
    }
    if (step == steps) break;

    if (!FLAGS_disable_wake) {
      auto shed = internal::ShedProcess(dt);
      internal::AdvectProcess(dt);
      internal::MorphingProcess(t);
      internal::AppendShedProcess(&shed);

      // TODO remove rings
      rings.wake_vortices().resize(vortices->size() - wake_offset);
      std::copy(vortices->begin() + wake_offset, vortices->end(),
                rings.wake_vortices().begin());
    }
  }
}

}  // namespace simulator
}  // namespace UVLM

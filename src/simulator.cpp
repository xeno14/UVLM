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
DEFINE_bool(use_joukowski, false, "force calculation method");

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
  if (morphings.size() == 0)
    LOG(FATAL) << FLAGS_run_name << " "
               << "No morphing";
  if (wings.size() == 0)
    LOG(FATAL) << FLAGS_run_name << " "
               << "No wing";

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

void AdvectProcess(std::vector<UVLM::VortexRing>* wake, const double dt) {
  const std::size_t wake_offset = WakeOffset();
  wake->clear();
  // copy current wake vortices
  wake->insert(wake->begin(), vortices->cbegin() + wake_offset,
               vortices->cend());
  UVLM::AdvectWake(wake, vortices->cbegin(), vortices->cend(), rings, inlet,
                   dt);
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
  auto gamma = ::UVLM::SolveLinearProblem(containers, morphings, inlet, t);
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
  for (std::size_t i = 0; i < containers.size(); i++) {
    const auto& morphing = morphings[i];
    for (const auto& v : containers[i]) {
      auto* mv = snapshot.add_morphing_velocities();
      mv->mutable_freestream()->CopyFrom(UVLM::Vector3dToPoint(inlet));
      Eigen::Vector3d vc;
      morphing.Velocity(&vc, v.ReferenceCentroid(), t);
      mv->mutable_center()->CopyFrom(::UVLM::Vector3dToPoint(vc));
      std::vector<Eigen::Vector3d> nodes;
      v.ForEachSegment([&](const auto& start, const auto& end,
                           const auto& start0, const auto& end0) {
        Eigen::Vector3d vel;
        morphing.Velocity(&vel, (start0 + end0) / 2, t);
        nodes.push_back(vel);
      });
      for (const auto& vn : nodes) {
        mv->add_nodes()->CopyFrom(UVLM::Vector3dToPoint(vn));
      }
    }
  }
  char filename[256];
  sprintf(filename, "%s/%08lu", output_path.c_str(), step);
  std::ofstream ofs(filename, std::ios::binary);
  CHECK((bool)ofs) << "Unable to open " << filename;
  snapshot.SerializeToOstream(&ofs);
}

// TODO multiple output
void CalcLoadProcess(const double t, const double dt) {
  std::vector<Eigen::Vector3d> loads;
  std::vector<double> data;
  data.push_back(t);
  if (FLAGS_use_joukowski) {
    LOG(INFO) << "joukowski";
  } else {
    LOG(INFO) << "Katz and Plotkin";
  }
  for (std::size_t i = 0; i < containers.size(); i++) {
    const auto& c = containers[i];
    const auto& c_prev = containers_prev[i];
    const auto& m = morphings[i];
    const double rho = 1;
    const double S = c.chord() * c.span();

    UVLM::calc_load::AerodynamicLoad load;
    if (FLAGS_use_joukowski) {
      load = UVLM::calc_load::CalcLoadJoukowski(c, c_prev, rings, m, inlet, rho,
                                                t, dt);
    } else {
      load = UVLM::calc_load::CalcLoadKatzPlotkin(c, c_prev, m, rings, inlet,
                                                  rho, t, dt);
      // auto load2  = UVLM::calc_load::CalcLoadJoukowski(c, c_prev, rings, m, inlet, rho,
      //                                           t, dt);
      // load.F.x() = load2.F.x();
    }
    const double U = inlet.norm();
    auto coeff = load.F / (0.5 * rho * U * U * S);

    data.push_back(coeff.x());
    data.push_back(coeff.y());
    data.push_back(coeff.z());
    data.push_back(load.Pin);
    data.push_back(load.Pout);
  }
  std::ofstream ofs(output_load_path, std::ios::app);
  CHECK((bool)ofs) << "Unable to open " << output_load_path;
  ofs << UVLM::util::join("\t", data.begin(), data.end()) << std::endl;
}

}  // namespace internal

void InitSimulator() {
  // do nothing
}

void AddWing(const ::UVLM::proto::Wing& wing, const ::UVLM::Morphing& m) {
  wings.emplace_back(wing);
  morphings.emplace_back(m);
  morphings.rbegin()->set_origin(UVLM::PointToVector3d(wing.origin()));
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
    LOG(INFO) << FLAGS_run_name << " "
              << "step: " << step;
    const double t = dt * step;
    const std::size_t wake_offset = internal::WakeOffset();

    // Save circulations of bound vortices at the previous step
    LOG(INFO) << "copy";
    containers_prev.resize(containers.size());
    CopyContainers(containers.begin(), containers.end(),
                   containers_prev.begin());

    // 連立方程式を解いて翼の上の循環を求める
    LOG(INFO) << "linear";
    internal::SolveLinearProblem(t);

    // TODO remove rings
    LOG(INFO) << "copy";
    CHECK(wake_offset == rings.bound_vortices().size())
        << "bound vortex size mismatch " << wake_offset << " "
        << rings.bound_vortices().size();
    std::copy(vortices->begin(), vortices->begin() + wake_offset,
              rings.bound_vortices().begin());
    internal::OutputSnapshot2(step, t);
    if (output_load_path.size()) {
      LOG(INFO) << "calc load";
      internal::CalcLoadProcess(t, dt);
    }
    if (step == steps) break;

    if (!FLAGS_disable_wake) {
      LOG(INFO) << "shed and advect";
      auto shed = internal::ShedProcess(dt);
      std::vector<UVLM::VortexRing> wake_next;
      internal::AdvectProcess(&wake_next, dt);
      internal::MorphingProcess(t);
      // Renew wake vortices
      LOG(INFO) << "copy";
      CHECK(wake_next.size() == vortices->size()-wake_offset) << "size mismatch";
      std::copy(wake_next.cbegin(), wake_next.cend(),
                vortices->begin() + wake_offset);
      internal::AppendShedProcess(&shed);

      // TODO remove rings
      rings.wake_vortices().resize(vortices->size() - wake_offset);
      CHECK((vortices->size() - wake_offset) == rings.wake_vortices().size())
          << "wake size mismatch";
      std::copy(vortices->begin() + wake_offset, vortices->end(),
                rings.wake_vortices().begin());
    }
  }
}

}  // namespace simulator
}  // namespace UVLM

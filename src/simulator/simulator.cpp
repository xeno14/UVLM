/**
 * @file simulator.cpp
 * @brief Add description here
 */

#include "simulator.h"

#include "../../proto/uvlm.pb.h"
#include "../calc_load/calc_load.h"
#include "../uvlm_vortex_ring.h"
#include "../morphing.h"
#include "../linear.h"
#include "../shed.h"
#include "../proto_adaptor.h"
#include "../util.h"
#include "../vortex_container.h"
#include "../wing_builder.h"

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

/**
 * @brief Connect trailing edge and newly shed wake
 */
void ConnectProcess() {
  // TODO multiple wings
  const std::size_t shed_size =
      std::distance(containers[0].edge_begin(), containers[0].edge_end());
  auto shed_first = vortices->begin() + WakeOffset();
  auto shed_last = shed_first + shed_size;
  auto edge = containers[0].edge_begin();

  if (shed_last > vortices->end()) return;

  for (auto it = shed_first; it != shed_last; ++it, ++edge) {
    it->nodes().resize(4);
    it->nodes()[0] = edge->nodes()[1];
    it->nodes()[3] = edge->nodes()[2];
  }
}

void MorphingProcess(const double t) {
  for (std::size_t i = 0; i < containers.size(); i++) {
    for (auto& vortex : containers[i]) {
      for (std::size_t j = 0; j < vortex.nodes().size(); j++) {
        morphings[i].Perfome(&vortex.nodes()[j], vortex.nodes0()[j], t);
      }
    }
  }
  ConnectProcess();
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
      load =
          UVLM::calc_load::CalcLoadJoukowski(c, c_prev, m, inlet, rho, t, dt);
    } else {
      load = UVLM::calc_load::CalcLoadKatzPlotkin(c, c_prev, m, rings, inlet,
                                                  rho, t, dt);
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

  // Main loop
  // See Katz and Plotkin p.408 Fig. 13.25
  // or Bueso2011 Fig 3.10
  for (std::size_t step = 0; step <= steps; ++step) {
    LOG(INFO) << FLAGS_run_name << " "
              << "step: " << step;
    const double t = dt * step;
    const std::size_t wake_offset = internal::WakeOffset();

    // Save circulations of bound vortices at the previous step
    LOG(INFO) << "copy";
    containers_prev.resize(containers.size());
    CopyContainers(containers.begin(), containers.end(),
                   containers_prev.begin());

    LOG(INFO) << "Kinematics";
    internal::MorphingProcess(t);

    LOG(INFO) << "Boundary Cond";
    internal::SolveLinearProblem(t);

    internal::OutputSnapshot2(step, t);
    if (output_load_path.size()) {
      LOG(INFO) << "Calc loads";
      internal::CalcLoadProcess(t, dt);
    }
    if (step == steps) break;

    if (!FLAGS_disable_wake) {
      LOG(INFO) << "Wake Rollup";
      // TODO multiple wings
      std::vector<UVLM::VortexRing> wake_next(containers[0].edge_begin(),
          containers[0].edge_end());
      wake_next.insert(wake_next.end(),
                       vortices->cbegin() + wake_offset, vortices->cend());
      UVLM::AdvectParallel(vortices->cbegin(), vortices->cend(),
                           wake_next.begin(), wake_next.end(),
                           inlet, dt);
      // update
      vortices->resize(wake_offset + wake_next.size());
      std::copy(wake_next.cbegin(), wake_next.cend(),
                vortices->begin() + wake_offset);
      // LOG(INFO) << vortices.size();
    }
  }
}

}  // namespace simulator
}  // namespace UVLM

/**
 * @file dump_centroid_main.cpp
 * @brief Add description here
 */

#include "../util.h"
#include "../../proto/uvlm.pb.h"
#include "../proto_adaptor.h"

#include <algorithm>
#include <fstream>
#include <gflags/gflags.h>
#include <glog/logging.h>
#include <glob.h>
#include <tuple>

DEFINE_string(input, "", "path to snapshot2");
DEFINE_string(output, "", "output path");

namespace {

std::string node_path;
std::string cent_path;
std::string vel_path;

inline auto ToVector(const Eigen::Vector3d& v) {
  return std::vector<double>{v.x(), v.y(), v.z()};
}

auto DumpCentroid(const UVLM::VortexContainer& c, std::ostream* os) {
  std::vector<Eigen::Vector3d> pos;
  for (const auto& v : c) {
    auto data = ToVector(v.Centroid());
    pos.push_back(v.Centroid());
    *os << UVLM::util::join("\t", data.begin(), data.end()) << "\n";
  }
  return pos;
}

void DumpNode(const UVLM::VortexContainer& c, std::ostream* os) {
  for (const auto& v : c) {
    std::vector<Eigen::Vector3d> nodes = v.nodes();
    nodes.push_back(nodes[0]);  // loop
    for (const auto& n : nodes) {
      auto data = ToVector(n);
      *os << UVLM::util::join("\t", data.begin(), data.end()) << "\n";
    }
    *os << "\n\n";
  }
}

void DumpVel(const std::vector<Eigen::Vector3d>& pos,
             const std::vector<Eigen::Vector3d>& vel, std::ostream* os) {
  if (pos.size() != vel.size()) return;
  for (std::size_t i = 0; i < pos.size(); i++) {
    auto data_pos = ToVector(pos[i]);
    auto data_vel = ToVector(vel[i]);
    *os << UVLM::util::join("\t", data_pos.begin(), data_pos.end()) << "\t" <<
      UVLM::util::join("\t", data_vel.begin(), data_vel.end()) << "\n";
  }
  *os << "\n\n";
}

std::string ResultName(const std::string& s) {
  auto pos = s.rfind("/");
  if (pos != std::string::npos) {
    return s.substr(pos + 1);
  }
  return s;
}

void Dump(const std::string& inpath, const std::string& outpath) {
  std::ifstream ifs(inpath, std::ios::binary);
  CHECK((bool)ifs) << "Unable to open " << FLAGS_input;

  std::ofstream ofs_node(node_path);
  std::ofstream ofs_cent(cent_path);
  std::ofstream ofs_vel(vel_path);
  CHECK(ofs_node) << "Unable to open";
  CHECK(ofs_cent) << "Unable to open";
  CHECK(ofs_vel) << "Unable to open";

  UVLM::proto::Snapshot2 snapshot;
  snapshot.ParseFromIstream(&ifs);

  std::vector<UVLM::VortexContainer> containers;
  auto vortices = UVLM::Snapshot2ToContainers(&containers, snapshot);

  std::vector<Eigen::Vector3d> vcenter, vfree;
  std::vector<std::vector<Eigen::Vector3d>> vnodes;
  UVLM::Snapshot2ToMorphingVelocities(&vcenter, &vnodes, &vfree, snapshot);

  std::vector<Eigen::Vector3d> pos;
  for (std::size_t index = 0; index < containers.size(); ++index) {
    DumpNode(containers[index], &ofs_node);
    auto postmp = DumpCentroid(containers[index], &ofs_cent);
    pos.insert(pos.end(), postmp.begin(), postmp.end());
  }
  DumpVel(pos, vcenter, &ofs_vel);
}
}  // anonymous namespace

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InstallFailureSignalHandler();
  FLAGS_logtostderr = true;

  auto res_name = ResultName(FLAGS_input);
  node_path = FLAGS_output + "/" + res_name + "_node.dat";
  cent_path = FLAGS_output + "/" + res_name + "_cent.dat";
  vel_path = FLAGS_output + "/" + res_name + "_vel.dat";

  Dump(FLAGS_input, FLAGS_output);

  printf(
      "splot '%s' w p, '%s' w l, '%s' usi 1:2:3:4:5:6 w vec\n"
      "pause -1\n",
      cent_path.c_str(), node_path.c_str(), vel_path.c_str());

  return 0;
}

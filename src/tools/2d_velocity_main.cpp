/**
 * @file 2d_velocity_main.cpp
 * @brief Add description here
 */

#include "../proto_adaptor.h"
#include "../shed.h"
#include "../util.h"

#include <cstdio>
#include <fstream>
#include <gflags/gflags.h>
#include <glog/logging.h>
#include <vector>
#include <map>
#include <set>

DEFINE_int32(resolution, 100, "number of points for each edge");
DEFINE_string(input, "", "snapshot2");
DEFINE_string(output, "", "output filename");
DEFINE_string(filetype, "csv", "csv or tsv");
DEFINE_string(plane, "xz", "coordinate plane. e.g. 'xy', 'yz', 'xz'");
DEFINE_string(range, "[-1:1]x[-1:1]@0",
              "Calculate range. e.g. if plane=='xy', [-1:1]x[-2:2]@0 will plot x "
              "in [-1, 1], y in [-2:2] at z=0");

namespace {

void CheckPlane(const std::string& plane) {
  CHECK(plane.size() == 2) << "FLAGS_plane must be size 2";
  CHECK(plane[0] != plane[1]) << "FLAGS_plane must have different charactors";
  CHECK(std::string("xyz").find(plane[0]) != std::string::npos)
      << "Charactors in FLAGS_plane must be one of 'xyz'";
  CHECK(std::string("xyz").find(plane[1]) != std::string::npos)
      << "Charactors in FLAGS_plane must be one of 'xyz'";
}

auto ParseRange(const std::string& range) {
  // TODO check format
  double x0, y0, x1, y1, cross_section;
  CHECK(sscanf(range.c_str(), "[%lf:%lf]x[%lf:%lf]@%lf", &x0, &y0, &x1, &y1,
               &cross_section) == 5)
      << "Something wrong with FLAGS_range format: " << range;
  return std::make_tuple(x0, y0, x1, y1, cross_section);
}

auto CreatePoints() {
  std::vector<Eigen::Vector3d> points;
  CheckPlane(FLAGS_plane);
  std::map<char, std::vector<double>> pos;
  std::set<char> keys {'x', 'y', 'z'};
  double min0, max0, min1, max1, cross_section;
  std::tie(min0, max0, min1, max1, cross_section) = ParseRange(FLAGS_range);
  LOG(INFO) << "range: " << min0 << "," << max0 << " x " << min1 << "," << max1
            << " at " << cross_section;

  // 0
  const char key0 = FLAGS_plane[0];
  pos[key0] =
      linspace(min0, max0, FLAGS_resolution);
  keys.erase(key0);
  // 1
  const char key1 = FLAGS_plane[1];
  pos[key1] =
      linspace(min1, max1, FLAGS_resolution);
  keys.erase(key1);
  // 2
  const char key2 = *keys.begin();
  pos[key2].push_back(cross_section);

  for (double x : pos['x']) {
    for (double y : pos['y']) {
      for (double z : pos['z']) {
        points.emplace_back(x, y, z);
      }
    }
  }
  return points;
}

auto GetVortices() {
  std::ifstream ifs(FLAGS_input, std::ios::binary);
  CHECK(ifs) << "Unable to open " << FLAGS_input;
  LOG(INFO) << "Input: " << FLAGS_input;

  UVLM::proto::Snapshot2 snapshot2;
  snapshot2.ParseFromIstream(&ifs);

  std::vector<UVLM::VortexContainer> containers;
  UVLM::Snapshot2ToContainers(&containers, snapshot2);

  CHECK(containers.size() > 0) << "no container";
  return containers.begin()->vortices();
}

struct Data {
  Eigen::Vector3d pos;
  Eigen::Vector3d vel;
};

template <class InputIterator1, class InputIterator2>
auto CalcData(InputIterator1 pos_first,
              InputIterator1 pos_last, InputIterator2 v_first,
              InputIterator2 v_last) {
  std::vector<Data> data;
  const std::size_t size = std::distance(pos_first, pos_last);
  data.resize(size);
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (std::size_t i = 0; i < size; i++) {
    const auto pos = *(pos_first + i);
    Eigen::Vector3d vel;
    UVLM::InducedVelocity(&vel, pos, v_first, v_last);
    data.at(i) = Data{pos, vel};
  }
  return data;
}

std::string Sep() {
  if (FLAGS_filetype == "csv") {
    return ",";
  } else if (FLAGS_filetype == "tsv") {
    return "\t";
  } else {
    LOG(FATAL) << "Invalid filetype: " << FLAGS_filetype;
  }
}


void Output(const std::vector<Data>& data, const std::string& output) {
  std::ofstream ofs(output);
  CHECK(ofs) << "Unable to open " << output;
  LOG(INFO) << "Output: " << output;
  const std::string sep = Sep();
  ofs << "#"
      << "x" << sep
      << "y" << sep
      << "z" << sep
      << "u" << sep
      << "v" << sep
      << "w" << std::endl;
  for (const auto& d : data) {
    ofs << d.pos.x() << sep
        << d.pos.y() << sep
        << d.pos.z() << sep
        << d.vel.x() << sep
        << d.vel.y() << sep
        << d.vel.z() << std::endl;
  }
}

void Write(const UVLM::proto::Snapshot2& snapshot, const std::string& output) {
  const auto points = CreatePoints();
  const auto vortices = GetVortices();
  const auto data = CalcData(points.cbegin(), points.cend(), vortices->cbegin(),
                             vortices->cend());
  Output(data, output);
}

}  // anonymous namespace

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InstallFailureSignalHandler();
  FLAGS_logtostderr = true;

  // std::vector<Eigen::Vector3d> points;
  // CreatePoints(&points);
  //
  // auto vortices = GetVortices();
  //
  // std::vector<Data> data;
  // CalcData(&data, points.cbegin(), points.cend(), vortices->cbegin(),
  //          vortices->cend());
  //
  // Output(data);

  return 0;
}

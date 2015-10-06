/**
 * @file 2d_velocity_main.cpp
 * @brief Add description here
 */

#include "../proto_adaptor.h"
#include "../shed.h"

#include <fstream>
#include <gflags/gflags.h>
#include <glog/logging.h>
#include <vector>

DEFINE_double(xmin, -1, "left most x");
DEFINE_double(zmin, -1.5, "down most z");
DEFINE_double(xmax, 3, "right most x");
DEFINE_double(zmax, 1.5, "up most z");
DEFINE_int32(resolution, 100, "number of points for each edge");
DEFINE_string(input, "", "snapshot2");
DEFINE_string(output, "", "output filename");
DEFINE_string(filetype, "csv", "csv or tsv");

namespace {

void CreatePoints(std::vector<Eigen::Vector3d>* points) {
  points->clear();
  const std::size_t resolution = FLAGS_resolution;
  const double xmin = FLAGS_xmin;
  const double zmin = FLAGS_zmin;
  const double xmax = FLAGS_xmax;
  const double zmax = FLAGS_zmax;
  const double dx = (xmax - xmin) / resolution;
  const double dz = (zmax - zmin) / resolution;

  for (std::size_t i = 0; i < resolution; i++) {
    for (std::size_t j = 0; j < resolution; j++) {
      double x = xmin + dx * i;
      double y = 0;
      double z = zmin + dz * j;
      points->emplace_back(x, y, z);
    }
  }
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
void CalcData(std::vector<Data>* data, InputIterator1 pos_first,
              InputIterator1 pos_last, InputIterator2 v_first,
              InputIterator2 v_last) {
  const std::size_t size = std::distance(pos_first, pos_last);
  data->resize(size);
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (std::size_t i = 0; i < size; i++) {
    const auto pos = *(pos_first + i);
    Eigen::Vector3d vel;
    UVLM::InducedVelocity(&vel, pos, v_first, v_last);
    data->at(i) = Data{pos, vel};
  }
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

void Output(const std::vector<Data>& data) {
  std::ofstream ofs(FLAGS_output);
  CHECK(ofs) << "Unable to open " << FLAGS_output;
  LOG(INFO) << "Output: " << FLAGS_output;
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

}  // anonymous namespace

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InstallFailureSignalHandler();
  FLAGS_logtostderr = true;

  std::vector<Eigen::Vector3d> points;
  CreatePoints(&points);

  auto vortices = GetVortices();

  std::vector<Data> data;
  CalcData(&data, points.cbegin(), points.cend(), vortices->cbegin(),
           vortices->cend());

  Output(data);

  return 0;
}

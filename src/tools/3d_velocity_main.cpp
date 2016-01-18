/**
 * @file 3d_velocity_main.cpp
 * @brief Add description here
 */
#include "../proto_adaptor.h"
#include "../recordio/recordio.h"
#include "../recordio/recordio_range.h"
#include "../shed.h"
#include "../util.h"

#include <boost/algorithm/string.hpp>
#include <cstdio>
#include <fstream>
#include <gflags/gflags.h>
#include <glob.h>
#include <glog/logging.h>
#include <map>
#include <set>
#include <string>
#include <vector>

DEFINE_int32(resolution, 100, "number of points for each edge");
DEFINE_string(input, "", "snapshot2 recordio");
DEFINE_string(output, "", "output filename");
DEFINE_string(range, "-1,-1,1,1", "min0,max0,min1,max1");
DEFINE_double(at, 0., "position of plane.");
DEFINE_string(prefix, "v", "prefix of output files");
DEFINE_int32(step, 0, "step");

namespace {

auto CreatePoints(const double xmin, const double xmax, const std::size_t Nx;
                  const double ymin, const double ymax, const std::size_t Ny;
                  const double zmin, const double zmax, const std::size_t Nz;) {
  const auto xs = linspace(xmin, xmax, Nx);
  const auto ys = linspace(ymin, ymax, Ny);
  const auto zs = linspace(zmin, zmax, Nz);

  for (double z : zs) {
    for (double y : ys) {
      for (double x : xs) {
        points.emplace_back(x, y, z);
      }
    }
  }
  return points;
}

std::vector<UVLM::VortexRing> GetVortices(const UVLM::proto::Snapshot2& snapshot) {
  std::vector<UVLM::VortexRing> res;
  res.resize(snapshot.vortices().size());
  std::transform(snapshot.vortices().begin(), snapshot.vortices().end(),
                 res.begin(),
                 [](const auto& v) { return UVLM::ProtoToVortexRing(v); });
  return res;
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


std::size_t Output(const std::vector<Data>& data, const std::string& output) {
  std::ofstream ofs(output);
  CHECK(ofs) << "Unable to open " << output;
  LOG(INFO) << "Output: " << output;
  const std::string sep = Sep();
  ofs << "x" << sep
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

std::size_t Write(const UVLM::proto::Snapshot2& snapshot, const std::string& output) {
  const auto points = CreatePoints();
  const auto vortices = GetVortices(snapshot);
  const auto data = CalcData(points.cbegin(), points.cend(), vortices.cbegin(),
                             vortices.cend());
  return Output(data, output);
}

std::string OutputPath(const std::string& output_path,
                       const std::size_t index) {
  char path[256];
  sprintf(path, "%s/%s%08lu.%s", output_path.c_str(), FLAGS_prefix.c_str(),
          index, FLAGS_filetype.c_str());
  return std::string(path);
}

std::size_t RecordioWriter(const std::string& input_path,
                    const std::string& output_path) {
  std::size_t index = 0;
  UVLM::proto::Snapshot2 snapshot;
  for (const auto& snapshot :
       recordio::ReaderRange<UVLM::proto::Snapshot2>(input_path)) {
    if (index == FLAGS_step) {
      return Write(snapshot, OutputPath(output_path, index));
    }
    ++index;
  }
  return 0;
}

}  // anonymous namespace

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InstallFailureSignalHandler();
  FLAGS_logtostderr = true;

  auto num = RecordioWriter(FLAGS_input, FLAGS_output);
  LOG(INFO) << "write " << num << " points";

  return 0;
}

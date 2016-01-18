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
DEFINE_double(xmin, 0, "min value of x");
DEFINE_double(xmax, 8, "max value of x");
DEFINE_double(ymin, -4, "min value of y");
DEFINE_double(ymax, 4, "max value of y");
DEFINE_double(zmin, -4, "min value of z");
DEFINE_double(zmax, 4, "max value of z");
DEFINE_string(prefix, "v", "prefix of output files");
DEFINE_uint64(step, 0, "step");

namespace {

auto CreatePoints(const double xmin, const double xmax, const std::size_t Nx,
                  const double ymin, const double ymax, const std::size_t Ny,
                  const double zmin, const double zmax, const std::size_t Nz) {
  const auto xs = linspace(xmin, xmax, Nx);
  const auto ys = linspace(ymin, ymax, Ny);
  const auto zs = linspace(zmin, zmax, Nz);

  std::vector<Eigen::Vector3d> res;
  for (double z : zs) {
    for (double y : ys) {
      for (double x : xs) {
        res.emplace_back(x, y, z);
      }
    }
  }
  return res;
}

std::vector<UVLM::VortexRing> GetVortices(
    const UVLM::proto::Snapshot2& snapshot) {
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

auto CalcData(const std::vector<Eigen::Vector3d>& points,
              const std::vector<UVLM::VortexRing>& vortices) {
  std::vector<Data> data(points.size());
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (std::size_t i = 0; i < points.size(); i++) {
    const auto& pos = points[i];
    Eigen::Vector3d vel;
    UVLM::InducedVelocity(&vel, pos, vortices.cbegin(), vortices.cend());
    vel += Eigen::Vector3d(1, 0, 0);  // freestream
    data.at(i) = Data{pos, vel};
  }
  return data;
}

std::size_t Write(const UVLM::proto::Snapshot2& snapshot,
                  const std::string& output) {
  const std::size_t Nx = FLAGS_resolution;
  const std::size_t Ny = FLAGS_resolution;
  const std::size_t Nz = FLAGS_resolution;
  const double xmin = FLAGS_xmin;
  const double xmax = FLAGS_xmax;
  const double ymin = FLAGS_ymin;
  const double ymax = FLAGS_ymax;
  const double zmin = FLAGS_zmin;
  const double zmax = FLAGS_zmax;
  const double dx = (xmax - xmin) / (Nx - 1);
  const double dy = (ymax - ymin) / (Ny - 1);
  const double dz = (zmax - zmin) / (Nz - 1);

  LOG(INFO) << "[" << xmin << "," << xmax << "]"
            << "x"
            << "[" << ymin << "," << ymax << "]"
            << "x"
            << "[" << zmin << "," << zmax << "]"
            << " dimensions "
            << "(" << Nx << "," << Ny << "," << Nz << ")";

  const auto points =
      CreatePoints(xmin, xmax, Nx, ymin, ymax, Ny, zmin, zmax, Nz);
  const auto vortices = GetVortices(snapshot);
  const auto result = CalcData(points, vortices);

  FILE* fp = fopen(output.c_str(), "w");
  CHECK(fp != NULL) << "Unable to open " << output;

  fprintf(fp, "# vtk DataFile Version 2.0\n");
  fprintf(fp, "HOGE FLOW\n");
  fprintf(fp, "ASCII\n");
  fprintf(fp, "DATASET STRUCTURED_POINTS\n");
  fprintf(fp, "DIMENSIONS %lu %lu %lu\n", Nx, Ny, Nz);
  fprintf(fp, "ORIGIN %f %f %f\n", xmin, ymin, zmin);
  fprintf(fp, "SPACING %f %f %f\n", dx, dy, dz);
  fprintf(fp, "POINT_DATA %lu\n", Nx * Ny * Nz);
  fprintf(fp, "VECTORS volume_vectors double\n");
  for (const auto& data : result) {
    const auto u = data.vel;
    fprintf(fp, "%lf %lf %lf\n", u.x(), u.y(), u.z());
  }
  fclose(fp);
  return result.size();
}

std::size_t RecordioWriter(const std::string& input_path,
                           const std::string& output_path) {
  std::size_t index = 0;
  UVLM::proto::Snapshot2 snapshot;
  for (const auto& snapshot :
       recordio::ReaderRange<UVLM::proto::Snapshot2>(input_path)) {
    if (index == FLAGS_step) {
      LOG(INFO) << "Step: " << FLAGS_step;
      return Write(snapshot, output_path);
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

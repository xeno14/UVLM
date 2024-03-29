/**
 * @file 2d_velocity_main.cpp
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
DEFINE_string(input, "", "snapshot2");
DEFINE_string(output, "", "output filename");
DEFINE_string(filetype, "csv", "csv or tsv");
DEFINE_string(plane, "xz", "coordinate plane. e.g. 'xy', 'yz', 'xz'");
DEFINE_string(range, "-1,-1,1,1", "min0,max0,min1,max1");
DEFINE_double(at, 0., "position of plane.");
DEFINE_bool(recordio, false, "use recordio instead of glob");
DEFINE_string(prefix, "v", "prefix of output files");

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
  std::array<double, 4> res;
  std::vector<std::string> splited;
  boost::algorithm::split(splited, FLAGS_range, boost::is_any_of(","));
  CHECK(splited.size() == 4);
  std::transform(splited.begin(), splited.end(), res.begin(),
                 [](const std::string& s) {
                   std::string::size_type sz;
                   return std::stod(s, &sz);
                 });
  return std::make_tuple(res[0], res[1], res[2], res[3]);
}

auto CreatePoints() {
  std::vector<Eigen::Vector3d> points;
  CheckPlane(FLAGS_plane);
  std::map<char, std::vector<double>> pos;
  std::set<char> keys {'x', 'y', 'z'};
  double min0, max0, min1, max1, cross_section;
  std::tie(min0, max0, min1, max1) = ParseRange(FLAGS_range);
  cross_section = FLAGS_at;
  LOG(INFO) << "range: [" << min0 << "," << max0 << "] x [" << min1 << "," << max1
            << "] at " << cross_section;

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

void Write(const UVLM::proto::Snapshot2& snapshot, const std::string& output) {
  const auto points = CreatePoints();
  const auto vortices = GetVortices(snapshot);
  const auto data = CalcData(points.cbegin(), points.cend(), vortices.cbegin(),
                             vortices.cend());
  Output(data, output);
}

std::string OutputPath(const std::string& output_path,
                       const std::size_t index) {
  char path[256];
  sprintf(path, "%s/%s%08lu.%s", output_path.c_str(), FLAGS_prefix.c_str(),
          index, FLAGS_filetype.c_str());
  return std::string(path);
}

std::size_t GlobWriter(const std::string& input_path,
                       const std::string& output_path) {
  glob_t globbuf;

  CHECK(glob(input_path.c_str(), 0, NULL, &globbuf) == 0);
  std::size_t index = 0;
  for (index = 0; index < globbuf.gl_pathc; index++) {
    std::ifstream in;
    CHECK((in.open(globbuf.gl_pathv[index], std::ios::binary), in));

    UVLM::proto::Snapshot2 snapshot;
    CHECK(snapshot.ParseFromIstream(&in));

    Write(snapshot, OutputPath(output_path, index));
  }
  globfree(&globbuf);
  return index;
}

std::size_t RecordioWriter(const std::string& input_path,
                    const std::string& output_path) {
  std::size_t index = 0;
  UVLM::proto::Snapshot2 snapshot;
  for (const auto& snapshot :
       recordio::ReaderRange<UVLM::proto::Snapshot2>(input_path)) {
    Write(snapshot, OutputPath(output_path, index));
    ++index;
  }
  return index;
}

}  // anonymous namespace

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InstallFailureSignalHandler();
  FLAGS_logtostderr = true;

  std::size_t num = 0;
  if (FLAGS_recordio) {
    num = RecordioWriter(FLAGS_input, FLAGS_output);
  } else {
    num = GlobWriter(FLAGS_input, FLAGS_output);
  }
  LOG(INFO) << "write " << num << " files";

  return 0;
}

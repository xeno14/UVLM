/**
 * @file velocity_at_main.cpp
 * @brief Add description here
 */

#include "../proto_adaptor.h"
#include "../shed.h"
#include "../util.h"
#include "../recordio/recordio_range.h"
#include "../csv_writer/csv_writer.h"

#include <boost/range/iterator_range.hpp>
#include <cstdio>
#include <fstream>
#include <gflags/gflags.h>
#include <glog/logging.h>
#include <vector>
#include <map>
#include <set>

DEFINE_string(input, "", "vortex sheet");
DEFINE_string(output, "", "output");
DEFINE_double(x, 0, "x");
DEFINE_double(y, 0, "y");
DEFINE_double(z, 0, "z");
DEFINE_int32(id, 0, "wing id.");
DEFINE_bool(wing, false, "include induced velocity due to wing");
DEFINE_bool(wake, false, "include induced velocity due to wake");

namespace {

inline auto GetWingSheet(const UVLM::proto::AllVortexSheets& sheet_proto) {
  return UVLM::proto_adaptor::FromVortexSheet(sheet_proto.wing());
}

inline auto GetWakeSheet(const UVLM::proto::AllVortexSheets& sheet_proto) {
  return UVLM::proto_adaptor::FromVortexSheet(sheet_proto.wake());
}

inline Eigen::Vector3d CalcVelocity(const Eigen::Vector3d& pos,
    const int id,
    const multiple_sheet::MultipleSheet<UVLM::VortexRing>& sheet) {
  Eigen::Vector3d res;
  UVLM::InducedVelocity(&res, pos, sheet.iterator_at(id, 0, 0),
                        sheet.iterator_at(id + 1, 0, 0));
  return res;
}

inline auto GetPos() { return Eigen::Vector3d(FLAGS_x, FLAGS_y, FLAGS_z); }

}  // anonymous namespace

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InstallFailureSignalHandler();
  FLAGS_logtostderr = true;

  CHECK(FLAGS_wing | FLAGS_wake)
      << "At least one of --wing or --wake must be true.";
  CHECK(FLAGS_id >= 0);

  std::ofstream ofs(FLAGS_output);
  CHECK(ofs);

  csv_writer::CSVWriter tsv(&ofs, "\t");
  tsv.write("step", "t", "x", "y", "z");

  const auto pos = GetPos();
  std::size_t count = 0;

  LOG(INFO) << "Velocity at " << pos.transpose() << " due to "
            << (FLAGS_wing ? "wing" : "")
            << (FLAGS_wing & FLAGS_wake ? " and " : "")
            << (FLAGS_wake ? "wake" : "")
            << " " << FLAGS_id;


  for (const auto& sheet_proto :
       recordio::ReaderRange<UVLM::proto::AllVortexSheets>(FLAGS_input)) {

    const Eigen::Vector3d v_wing =
        !FLAGS_wing ? Eigen::Vector3d::Zero()
                    : CalcVelocity(pos, FLAGS_id, GetWingSheet(sheet_proto));
    const Eigen::Vector3d v_wake =
        !FLAGS_wake ? Eigen::Vector3d::Zero()
                    : CalcVelocity(pos, FLAGS_id, GetWakeSheet(sheet_proto));
    const Eigen::Vector3d result = v_wing + v_wake;
    const double t = sheet_proto.t();

    tsv.write(count, t, result.x(), result.y(), result.z());
    ++count;
  }
  tsv.flush();
  LOG(INFO) << count << " lines are written."
            << (count == 0 ? " Is input file correct?" : "");

  return 0;
}

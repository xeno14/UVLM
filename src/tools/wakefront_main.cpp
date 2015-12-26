/**
 * @file wakefront_main.cpp
 * @brief Add description here
 */

#include "../util.h"
#include "../recordio/recordio.h"
#include "../proto_adaptor.h"

#include <iostream>
#include <cstdio>
#include <fstream>
#include <gflags/gflags.h>
#include <glog/logging.h>
#include <vector>
#include <map>
#include <set>

DEFINE_string(input, "", "AllVortexSheets recordio");
DEFINE_string(output, "", "output filename");
DEFINE_int32(start, 0, "start tracking wake front");
DEFINE_int32(id, 0, "wing id");

namespace {

int WakeFront() {
  std::ifstream ifs(FLAGS_input, std::ios::binary);
  CHECK(ifs) << "Unable to open " << FLAGS_input;
  LOG(INFO) << "Input: " << FLAGS_input;
  recordio::RecordReader reader(&ifs);

  std::ofstream ofs(FLAGS_output);
  CHECK(ofs) << "Unable to open " << FLAGS_output;
  LOG(INFO) << "Output: " << FLAGS_output;

  std::size_t steps = 0;
  std::size_t count = 0;
  UVLM::proto::AllVortexSheets sheets;
  while (reader.ReadProtocolMessage(&sheets)) {
    const double t = sheets.t();
    const auto wake = UVLM::proto_adaptor::FromVortexSheet(sheets.wake());

    const long target_row = wake.rows() - FLAGS_start - 1;

    if (target_row < 0) continue;

    ofs << t;

    const auto first = wake.iterator_at(FLAGS_id, target_row, 0);
    const auto last = wake.iterator_at(FLAGS_id, target_row + 1, 0);
    for (auto it = first; it!=last; ++it) {
      ofs << "\t" << it->nodes().begin()->x();
    }
    ofs << "\t" << (last - 1)->nodes()[1].x();
    ofs << std::endl;
    
    ++steps;
    ++count;
  }
  LOG(INFO) << count << " lines are saved";

  return 0;
}

}  // anonymous


int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InstallFailureSignalHandler();
  FLAGS_logtostderr = true;

  return WakeFront();
}

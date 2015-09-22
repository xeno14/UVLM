/**
 * @file calc_load.cpp
 * @brief Add description here
 */

#include "calc_load.h"

#include <fstream>
#include <glob.h>
#include <glog/logging.h>

namespace UVLM {
namespace calc_load {
namespace internal {

proto::Snapshot2 ReadSnapshot2(const std::string& filename) {
  std::ifstream ifs(filename, std::ios::binary);
  CHECK((bool)ifs) << "Unable to open " << filename;
  proto::Snapshot2 snapshot2;
  snapshot2.ParseFromIstream(&ifs);
  return snapshot2;
}

std::vector<proto::Snapshot2> ReadSnapshot2(
    const std::vector<std::string>& filenames) {
  std::vector<proto::Snapshot2> res;
  for (const auto& filename : filenames) {
    res.push_back(ReadSnapshot2(filename.c_str()));
  }
  std::sort(res.begin(), res.end(),
            [](const auto& l, const auto& r) { return l.t() < r.t(); });
  return res;
}

}  // namespace internal

std::vector<std::string> Snapshot2Paths(const std::string& pattern) {
  glob_t globbuf;

  int ret = glob(pattern.c_str(), 0, NULL, &globbuf);
  if (ret == GLOB_NOMATCH) {
    LOG(FATAL) << "Pattern no match";
  }

  std::vector<std::string> res;
  for (std::size_t i = 0; i < globbuf.gl_pathc; i++) {
    res.emplace_back(globbuf.gl_pathv[i]);
  }
  globfree(&globbuf);

  std::sort(res.begin(), res.end());
  return res;
}

void Start(const ::UVLM::Morphing& m, const std::string& output_path) {
  auto filepaths = Snapshot2Paths(output_path + "/*");

  if (filepaths.size() < 2) return;

}

}  // namespace calc_load
}  // namespace UVLM

/**
 * @file calc_load.cpp
 * @brief Add description here
 */

#include "../proto_adaptor.h"
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

/**
 * @param s0 snapshot2 at t-Δt
 * @param s1 snapshot2 at t
 */
std::vector<Eigen::Vector3d> Calc(const proto::Snapshot2& s0, const proto::Snapshot2& s1) {
  std::vector<UVLM::VortexContainer> c0, c1;
  UVLM::Snapshot2ToContainers(&c0, s0); UVLM::Snapshot2ToContainers(&c1, s1);

  if (c0.size() != c1.size()) {
    LOG(FATAL) << "Container size not match";
  }

  const double t = s1.t();
  const double dt = s1.t() - s0.t();
  const double RHO = 1.0;
  const double U = 1;
  const double ALPHA = 0 * M_PI / 180;
  UVLM::Morphing morphing;  // do nothing...
  Eigen::Vector3d inlet(U * cos(ALPHA), 0, U * sin(ALPHA)); // TODO use inlet inside of morphing

  std::vector<Eigen::Vector3d> res;

  auto wake_iterator = GetWake(c1);
  for (std::size_t i = 0; i < c1.size(); i++) {
    auto load = CalcLoad(c1[i], c0[i], wake_iterator.first,
                         wake_iterator.second, morphing, inlet, RHO, t, dt);
    load /= (0.5 * RHO * U * U * c1[i].chord() * c1[i].span());
    res.emplace_back(load);
  }
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

/**
 * @todo morphing should be list.
 */
void Start(const ::UVLM::Morphing& m,
           const std::string& output_path) {
  auto filepaths = Snapshot2Paths(output_path + "/*");

  if (filepaths.size() < 2) return;
  // snapshot at t-dt
  proto::Snapshot2 s0 = internal::ReadSnapshot2(*(filepaths.rbegin() + 1));
  // snapshot at t
  proto::Snapshot2 s1 = internal::ReadSnapshot2(*filepaths.rbegin());

  auto result = internal::Calc(s0, s1);
}

}  // namespace calc_load
}  // namespace UVLM

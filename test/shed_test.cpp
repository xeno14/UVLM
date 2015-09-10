#include <gtest/gtest.h>

#include "shed.h"
#include "testutil.h"

using UVLM::UVLMVortexRing;
using UVLM::VortexRing;


TEST(InducedVelocityTest, many_rings) {
  VortexRing v;
  v.set_gamma(4 * M_PI);
  v.PushNode(1, 1, 0).PushNode(-1, 1, 0).PushNode(-1, -1, 0).PushNode(1, -1, 0);

  std::vector<VortexRing> rings;
  const std::size_t N = 20;
  for (std::size_t i = 0; i < N; i++) {
    rings.emplace_back(v);
  }

  Eigen::Vector3d pos(0, 0, 0);

  // 4sqrt(2)の流れがN個重ねあわせる
  Eigen::Vector3d res;
  InducedVelocity(&res, pos, rings.begin(), rings.end());
  EXPECT_VECTOR3D_EQ(0, 0, N * 4 * M_SQRT2, res);
}

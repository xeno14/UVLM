#include <gtest/gtest.h>

#include "shed.h"
#include "testutil.h"

using UVLM::UVLMVortexRing;
using UVLM::VortexRing;

class InducedVelocityTest : public ::testing::Test {
 protected:
  InducedVelocityTest() {
    v.set_gamma(4 * M_PI);
    v.PushNode(1, 1, 0)
      .PushNode(-1, 1, 0)
      .PushNode(-1, -1, 0)
      .PushNode(1, -1, 0);
  }
  VortexRing v;
};

TEST_F(InducedVelocityTest, many_rings) {
  std::vector<VortexRing> rings;
  const std::size_t N = 20;
  for (std::size_t i = 0; i < N; i++) {
    rings.emplace_back(v);
  }
  ASSERT_EQ(20, rings.size());

  Eigen::Vector3d pos(0, 0, 0);

  // -4sqrt(2)の流れがN個重ねあわせる
  Eigen::Vector3d res;
  InducedVelocity(&res, pos, rings.begin(), rings.end());
  EXPECT_VECTOR3D_EQ(0, 0, -4 * M_SQRT2 * N, res);
}

TEST_F(InducedVelocityTest, chordwise_many_rings) {
  std::vector<VortexRing> rings;
  const std::size_t N = 20;
  for (std::size_t i = 0; i < N; i++) {
    rings.emplace_back(v);
  }

  Eigen::Vector3d pos(0, 0, 0);

  // 2sqrt(2)の流れがN個重ねあわせる
  Eigen::Vector3d res;
  ChordwiseInducedVelocity(&res, pos, rings.begin(), rings.end());
  EXPECT_VECTOR3D_EQ(0, 0, -4 * M_SQRT2 / 2 * N, res);
}

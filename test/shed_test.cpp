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

TEST(ConnectTrailingEdge, test) {
  std::vector<VortexRing> edge{GetSquareRing(2, 1, 1), GetSquareRing(2, 1, 3)};
  std::vector<VortexRing> wake{GetSquareRing(2, 99, 1),
                               GetSquareRing(2, 99, 3)};
  UVLM::ConnectTrailingEdge(edge.cbegin(), edge.cend(), wake.begin());
  EXPECT_VECTOR3D_EQ(2, 0, 0, wake[0].nodes()[0]);
  EXPECT_VECTOR3D_EQ(100, 0, 0, wake[0].nodes()[1]);
  EXPECT_VECTOR3D_EQ(100, 2, 0, wake[0].nodes()[2]);
  EXPECT_VECTOR3D_EQ(2, 2, 0, wake[0].nodes()[3]);
  EXPECT_VECTOR3D_EQ(2, 2, 0, wake[1].nodes()[0]);
  EXPECT_VECTOR3D_EQ(100, 2, 0, wake[1].nodes()[1]);
  EXPECT_VECTOR3D_EQ(100, 4, 0, wake[1].nodes()[2]);
  EXPECT_VECTOR3D_EQ(2, 4, 0, wake[1].nodes()[3]);
}

TEST(Advect, no_vortex) {
  std::vector<VortexRing> wake{GetSquareRing(2, 0, 0)};
  std::vector<VortexRing> dummy;
  Eigen::Vector3d freestream(1, 0, 0);
  UVLM::Advect(dummy.begin(), dummy.end(), wake.begin(), wake.end(), freestream,
               0.1);
  auto nodes = wake.begin()->nodes();
  EXPECT_VECTOR3D_EQ(-1 + 0.1, -1, 0, nodes[0]);
  EXPECT_VECTOR3D_EQ( 1 + 0.1, -1, 0, nodes[1]);
  EXPECT_VECTOR3D_EQ( 1 + 0.1,  1, 0, nodes[2]);
  EXPECT_VECTOR3D_EQ(-1 + 0.1,  1, 0, nodes[3]);
}

TEST(Advect, parallel) {
  std::vector<VortexRing> wake{GetSquareRing(2, 0, 0)};
  std::vector<VortexRing> dummy;
  Eigen::Vector3d freestream(1, 0, 0);
  UVLM::AdvectParallel(dummy.begin(), dummy.end(), wake.begin(), wake.end(),
                       freestream, 0.1);
  auto nodes = wake.begin()->nodes();
  EXPECT_VECTOR3D_EQ(-1 + 0.1, -1, 0, nodes[0]);
  EXPECT_VECTOR3D_EQ( 1 + 0.1, -1, 0, nodes[1]);
  EXPECT_VECTOR3D_EQ( 1 + 0.1,  1, 0, nodes[2]);
  EXPECT_VECTOR3D_EQ(-1 + 0.1,  1, 0, nodes[3]);
}

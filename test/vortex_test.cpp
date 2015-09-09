#include <gtest/gtest.h>

#include "vortex.h"

const double EPS = 1e-10;

class DistanceLineAndPointTest : public ::testing::Test {
 protected:
  Eigen::Vector3d a, b, c;
};

TEST_F(DistanceLineAndPointTest, test1) {
  a << 0, 0, -1;
  b << 0, 0, 1;
  c << 1, 0, 0;
  EXPECT_NEAR(1, UVLM::DistanceLineAndPoint(a, b, c), EPS);
}

TEST_F(DistanceLineAndPointTest, test2) {
  a << 0, 0, 0;
  b << 1, 1, 0;
  c << 1, 0, 0;
  EXPECT_NEAR(M_SQRT2/2, UVLM::DistanceLineAndPoint(a, b, c), EPS);
}

class BiotSavartLawTest : public ::testing::Test {
  protected:
  Eigen::Vector3d a, b, c, result;
};

TEST_F(BiotSavartLawTest, infinite_length) {
  // When the line segment is sufficiently long, strength equals to 1/4pi
  a << 0, 0, -100000;
  b << 0, 0,  100000;
  c << 1, 0, 0;
  UVLM::BiotSavartLaw(&result, a, b, c);
  EXPECT_NEAR(0, result.x(), EPS);
  EXPECT_NEAR(1./2./M_PI, result.y(), EPS);
  EXPECT_NEAR(0, result.z(), EPS);
}

TEST_F(BiotSavartLawTest, finite_length) {
  a << 0, 0, -1;
  b << 0, 0, 1;
  c << 1, 0, 0;
  UVLM::BiotSavartLaw(&result, a, b, c);
  EXPECT_NEAR(0, result.x(), EPS);
  EXPECT_NEAR(M_SQRT2/4./M_PI, result.y(), EPS);
  EXPECT_NEAR(0, result.z(), EPS);
}


class VortexRingTest : public ::testing::Test {
 protected:
  VortexRingTest() {}
  ~VortexRingTest() {}
  virtual void SetUp() {}
  virtual void TearDown() { ring.Clear(); }
  UVLM::VortexRing ring;
};

TEST_F(VortexRingTest, Assemble) {
  ring.PushNode(Vector3d(0, 0, 0))
      .PushNode(Vector3d(1, 0, 0))
      .PushNode(Vector3d(1, 1, 0))
      .PushNode(Vector3d(0, 1, 0));
  const auto& nodes = ring.nodes();
  ASSERT_EQ(4, nodes.size());
  EXPECT_EQ(Vector3d(0, 0, 0), nodes[0]);
  EXPECT_EQ(Vector3d(1, 0, 0), nodes[1]);
  EXPECT_EQ(Vector3d(1, 1, 0), nodes[2]);
  EXPECT_EQ(Vector3d(0, 1, 0), nodes[3]);
}

TEST_F(VortexRingTest, BiotSavartLaw) {
  ring.set_gamma(4*M_PI);
  ring.PushNode(1, 1, 0)
      .PushNode(-1, 1, 0)
      .PushNode(-1, -1, 0)
      .PushNode(1, -1, 0);
  Vector3d result;
  Vector3d pos(0, 0, 0);
  ring.BiotSavartLaw(&result, pos);
  EXPECT_DOUBLE_EQ(0, result.x());
  EXPECT_DOUBLE_EQ(0, result.y());
  EXPECT_DOUBLE_EQ(4*M_SQRT2, result.z());
}
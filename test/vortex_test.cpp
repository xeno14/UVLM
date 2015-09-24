#include <gtest/gtest.h>

#include "vortex.h"
#include "testutil.h"

using UVLM::VortexRing;

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

  UVLM::VortexRing GetSquaredRing(double l) {
    UVLM::VortexRing res;
    res.PushNode(Vector3d(0, 0, 0))
       .PushNode(Vector3d(l, 0, 0))
       .PushNode(Vector3d(l, l, 0))
       .PushNode(Vector3d(0, l, 0));
    return res;
  }
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

TEST_F(VortexRingTest, ChordwiseBiotSavartLaw) {
  ring.set_gamma(4*M_PI);
  ring.PushNode(1, 1, 0)
      .PushNode(-1, 1, 0)
      .PushNode(-1, -1, 0)
      .PushNode(1, -1, 0);
  Vector3d result;
  Vector3d pos(0, 0, 0);
  ring.ChordwiseBiotSavartLaw(&result, pos);
  EXPECT_DOUBLE_EQ(0, result.x());
  EXPECT_DOUBLE_EQ(0, result.y());
  EXPECT_DOUBLE_EQ(2*M_SQRT2, result.z());
}

TEST_F(VortexRingTest, Normal) {
  auto v = GetSquaredRing(20000);
  auto n = v.Normal();
  EXPECT_VECTOR3D_EQ(0, 0, 1, n);
}

TEST_F(VortexRingTest, Tangent) {
  auto v = GetSquaredRing(20000);
  auto t = v.Tangent();
  EXPECT_VECTOR3D_EQ(1, 0, 0, t);
}

TEST_F(VortexRingTest, AngleOfAttack) {
  auto v = GetSquaredRing(1);
  const double alpha = 0.1;
  Vector3d Q = {cos(alpha), 0 , sin(alpha)};
  EXPECT_DOUBLE_EQ(alpha, v.AngleOfAttack(Q));
}


class ChordSpanTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
    v.Clear();
    v.PushNode({0, 0, 0});
    v.PushNode({1, 0, 0});
    v.PushNode({1, 2.1, 0});
    v.PushNode({0, 2, 0});
  }

  VortexRing v;
};

TEST_F(ChordSpanTest, tanvec_c) {
  EXPECT_VECTOR3D_EQ(1, 0, 0, v.TanVecChord());
}

TEST_F(ChordSpanTest, tanvec_b) {
  EXPECT_VECTOR3D_EQ(0, 1, 0, v.TanVecSpan());
}

TEST_F(ChordSpanTest, c) {
  EXPECT_DOUBLE_EQ(1, v.CalcC());
}

TEST_F(ChordSpanTest, b) {
  EXPECT_DOUBLE_EQ(2, v.CalcB());
}

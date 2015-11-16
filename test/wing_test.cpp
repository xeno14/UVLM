#include <gtest/gtest.h>

#include "testutil.h"
#include "wing/wing.h"

using UVLM::wing::RectGenerator;
using UVLM::wing::WholeWing;
using UVLM::proto::Wing;
using UVLM::PointsToVector;

namespace {
const double EPS=1e-10;
}

/*
 * wing: AR=4 (chord=1, span=4), rows=2, cols=4
 * half: chord=1, span=2, rows=2, cols=2
 */
class WholeWingTest : public ::testing::Test {
 protected:
  WholeWingTest() : rect(1, 2, 2, 2) {}
  virtual void SetUp() {
    wing_base.Clear();
    rect.Generate(&wing_base);
  }
  RectGenerator rect;
  Wing wing_base;
};

TEST_F(WholeWingTest, mirror) {
  Wing wing;
  WholeWing(&wing, wing_base);
  auto pos = PointsToVector(wing.points());
  ASSERT_EQ(15, pos.size());
  EXPECT_VECTOR3D_NEAR(0, -2, 0, pos[0], EPS);
  EXPECT_VECTOR3D_NEAR(0, -1, 0, pos[1], EPS);
  EXPECT_VECTOR3D_NEAR(0, 0, 0, pos[2], EPS);
  EXPECT_VECTOR3D_NEAR(0, 1, 0, pos[3], EPS);
  EXPECT_VECTOR3D_NEAR(0, 2, 0, pos[4], EPS);
  EXPECT_VECTOR3D_NEAR(0.5, -2, 0, pos[5], EPS);
  EXPECT_VECTOR3D_NEAR(1, -2, 0, pos[10], EPS);
  EXPECT_VECTOR3D_NEAR(1, 2, 0, pos[14], EPS);
}

TEST_F(WholeWingTest, origin) {
  UVLM::wing::SetOrigin(&wing_base, {1, 2, 3});
  Wing wing;
  WholeWing(&wing, wing_base);
  auto pos = PointsToVector(wing.points());
  EXPECT_VECTOR3D_NEAR(1, 0, 3, pos[0], EPS);
  EXPECT_VECTOR3D_NEAR(1, 1, 3, pos[1], EPS);
  EXPECT_VECTOR3D_NEAR(1, 2, 3, pos[2], EPS);
  EXPECT_VECTOR3D_NEAR(1, 3, 3, pos[3], EPS);
  EXPECT_VECTOR3D_NEAR(1, 4, 3, pos[4], EPS);
  EXPECT_VECTOR3D_NEAR(1.5, 0, 3, pos[5], EPS);
  EXPECT_VECTOR3D_NEAR(2, 0, 3, pos[10], EPS);
  EXPECT_VECTOR3D_NEAR(2, 4, 3, pos[14], EPS);
}

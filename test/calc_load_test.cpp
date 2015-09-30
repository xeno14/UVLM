#include <gtest/gtest.h>

#include "../src/calc_load/calc_load.h"
#include "testutil.h"

using UVLM::VortexRing;
using namespace UVLM::calc_load;

class CalcLoadFuncTest : public ::testing::Test {
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

TEST_F(CalcLoadFuncTest, dummy) {

}

class LocalUnitVectorTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
    n << 0, 0, 1;
    t << 1, 0, 0;
  }
  Eigen::Vector3d n, t;
  Eigen::Vector3d e_lift, e_drag;
};

TEST(ProjectionOperatorTest, CalcProjectionOperator) {
  Eigen::Vector3d Um(2, 0, 0);
  Eigen::Matrix3d P = internal::CalcProjectionOperator(Um);
  Eigen::Vector3d result;
  result = P.block(0, 0, 3, 1);
  EXPECT_VECTOR3D_EQ(0, 0, 0, result);
  result = P.block(0, 1, 3, 1);
  EXPECT_VECTOR3D_EQ(0, 1, 0, result);
  result = P.block(0, 2, 3, 1);
  EXPECT_VECTOR3D_EQ(0, 0, 1, result);
}

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

TEST_F(LocalUnitVectorTest, MatrixForLocalUnitVector) {
  Eigen::Matrix3d m;
  UVLM::MatrixForLocalUnitVector(&m, n, t, 0);
  EXPECT_EQ(Eigen::Matrix3d::Identity(), m);
}

TEST_F(LocalUnitVectorTest, LocalUnitVector) {
  UVLM::LocalUnitVector(&e_lift, &e_drag, n, t, 0);
  EXPECT_VECTOR3D_EQ(0, 0, 1, e_lift);
  EXPECT_VECTOR3D_EQ(1, 0, 0, e_drag);
}

TEST_F(LocalUnitVectorTest, LocalUnitVector2) {
  const double alpha = M_PI / 6;
  UVLM::LocalUnitVector(&e_lift, &e_drag, n, t, alpha);
  EXPECT_VECTOR3D_EQ(-sin(alpha), 0, cos(alpha), e_lift);
  EXPECT_VECTOR3D_EQ(cos(alpha), 0, sin(alpha), e_drag);
}


TEST(ProjectionOperatorTest, GetProjectionOperator) {
  Eigen::Vector3d Um(2, 0, 0);
  Eigen::Matrix3d P = internal::GetProjectionOperator(Um);
  Eigen::Vector3d result;
  result = P.block(0, 0, 3, 1);
  EXPECT_VECTOR3D_EQ(0, 0, 0, result);
  result = P.block(0, 1, 3, 1);
  EXPECT_VECTOR3D_EQ(0, 1, 0, result);
  result = P.block(0, 2, 3, 1);
  EXPECT_VECTOR3D_EQ(0, 0, 1, result);
}

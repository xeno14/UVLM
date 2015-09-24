#include <gtest/gtest.h>

#include "../src/calc_load/calc_load.h"
#include "testutil.h"

using UVLM::VortexRing;

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
};

TEST_F(LocalUnitVectorTest, MatrixForLocalUnitVector) {
  Eigen::Matrix3d m;
  UVLM::MatrixForLocalUnitVector(&m, n, t, 0);
  EXPECT_EQ(Eigen::Matrix3d::Identity(), m);
}

TEST_F(LocalUnitVectorTest, LocalUnitVector) {
  Eigen::Vector3d e_lift, e_drag;
  UVLM::LocalUnitVector(&e_lift, &e_drag, n, t, 0);
  EXPECT_VECTOR3D_EQ(0, 0, 1, e_lift);
  EXPECT_VECTOR3D_EQ(1, 0, 0, e_drag);
}

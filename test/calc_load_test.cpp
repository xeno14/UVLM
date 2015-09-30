#include <gtest/gtest.h>

#include "../src/calc_load/calc_load.h"
#include "testutil.h"

using UVLM::VortexRing;
using namespace UVLM::calc_load;

namespace {
const double EPS = 1e-10;
}

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

TEST_F(CalcLoadFuncTest, CalcProjectionOperator) {
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

TEST_F(CalcLoadFuncTest, CalcUm) {
  UVLM::Morphing m;
  m.set_plug([](double t) { return sin(t); });
  Eigen::Vector3d freestream(1, 0, 0);
  EXPECT_VECTOR3D_NEAR(1, 0, -1, 
                       internal::CalcUm(m, {0, 0, 0}, freestream, 0), EPS);
  EXPECT_VECTOR3D_NEAR(1, 0, -0.5,
                       internal::CalcUm(m, {0, 0, 0}, freestream, M_PI / 3), EPS);
}

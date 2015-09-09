#include <gtest/gtest.h>

#include "calc_load.h"

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

TEST_F(CalcLoadFuncTest, c) {
  EXPECT_DOUBLE_EQ(1, UVLM::calc_load::CalcC(v));
}

TEST_F(CalcLoadFuncTest, b) {
  EXPECT_DOUBLE_EQ(2, UVLM::calc_load::CalcB(v));
}

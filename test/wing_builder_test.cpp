#include <gtest/gtest.h>

#include "wing_builder.h"

#include <iostream>

using namespace UVLM;

class WingBuilderTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
    vortices = std::make_shared<std::vector<VortexRing>>();
  }

  virtual void TearDown() {
    
  }
  std::shared_ptr<std::vector<VortexRing>> vortices;
};

TEST_F(WingBuilderTest, transform_row) {
  std::vector<Eigen::Vector3d> row{{1, 0, 2}, {1, 1, 2}, {1, 2, 2}};
  auto result = internal::TransfromRow(row.begin(), row.end());

  // 符号の反転の順番のチェック
  EXPECT_DOUBLE_EQ(-2, result[0].y());
  EXPECT_DOUBLE_EQ(-1, result[1].y());
  EXPECT_DOUBLE_EQ( 0, result[2].y());
  EXPECT_DOUBLE_EQ( 0, result[3].y());
  EXPECT_DOUBLE_EQ( 1, result[4].y());
  EXPECT_DOUBLE_EQ( 2, result[5].y());

  // x, zが不変なことをチェック
  for (const auto& r : result) {
    EXPECT_DOUBLE_EQ(1, r.x());
    EXPECT_DOUBLE_EQ(2, r.z());
  }
}

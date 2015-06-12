#include <gtest/gtest.h>

#include "wing_builder.h"

#include <iostream>

using namespace UVLM;

class WingBuilderTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
    vortices = std::make_shared<std::vector<VortexRing>>();
    containers.clear();

    int rows = 2, cols = 3;
    wing.set_rows(rows);
    wing.set_cols(cols);
    for (int i=0; i<=rows; i++) {
      for (int j=0; j<=cols; j++) {
        auto* p = wing.add_points();
        p->set_x(i); p->set_y(j); p->set_z(0);
      }
    }
  }

  std::shared_ptr<std::vector<VortexRing>> vortices;
  std::vector<VortexContainer> containers;
  proto::Wing wing;
};

TEST_F(WingBuilderTest, transform_row) {
  std::vector<Eigen::Vector3d> row{{1, 0, 2}, {1, 1, 2}, {1, 2, 2}};
  auto result = internal::TransfromRow(row.begin(), row.end());

  ASSERT_EQ(5, result.size());

  // 符号の反転の順番のチェック
  EXPECT_DOUBLE_EQ(-2, result[0].y());
  EXPECT_DOUBLE_EQ(-1, result[1].y());
  EXPECT_DOUBLE_EQ( 0, result[2].y());
  EXPECT_DOUBLE_EQ( 1, result[3].y());
  EXPECT_DOUBLE_EQ( 2, result[4].y());

  // x, zが不変なことをチェック
  for (const auto& r : result) {
    EXPECT_DOUBLE_EQ(1, r.x());
    EXPECT_DOUBLE_EQ(2, r.z());
  }
}

TEST_F(WingBuilderTest, total_size) {
  std::vector<internal::WingHolder> holders(3);
  holders[0].cols = holders[0].rows = 1;
  holders[1].cols = holders[1].rows = 2;
  holders[2].cols = holders[2].rows = 3;
  EXPECT_EQ(1+4+9, WingBuilder::CountTotalSize(holders));
}

TEST_F(WingBuilderTest, AddWing) {
  WingBuilder builder(&containers, vortices);
  const auto& holders = builder.holders();
  builder.AddWing(wing);
  
  ASSERT_EQ(1, holders.size());
  EXPECT_EQ(2, holders[0].rows);
  EXPECT_EQ(6, holders[0].cols);

  const auto& holder = holders[0];
  EXPECT_DOUBLE_EQ(0,  holder.points[0].x());
  EXPECT_DOUBLE_EQ(-3, holder.points[0].y());
  EXPECT_DOUBLE_EQ(0,  holder.points[0].z());
}

TEST_F(WingBuilderTest, add_multiple_wings) {
  WingBuilder builder(&containers, vortices);

  proto::Wing wing2;
  wing2.CopyFrom(wing);
  auto* origin = wing2.mutable_origin();
  origin->set_x(0); origin->set_y(0); origin->set_z(1);
  EXPECT_DOUBLE_EQ(1, wing2.origin().z());

  builder.AddWing(wing).AddWing(wing2);

  const auto& holders = builder.holders();
  ASSERT_EQ(2, holders.size());
  EXPECT_EQ(2, holders[1].rows);
  EXPECT_EQ(6, holders[1].cols);

  const auto& holder = holders[1];
  EXPECT_DOUBLE_EQ(0,  holder.points[0].x());
  EXPECT_DOUBLE_EQ(-3, holder.points[0].y());
  EXPECT_DOUBLE_EQ(1,  holder.points[0].z());
  EXPECT_DOUBLE_EQ(1, holder.points[0].z());
  EXPECT_DOUBLE_EQ(1, holder.points[1].z());
  EXPECT_DOUBLE_EQ(1, holder.points[6].z());
}

TEST_F(WingBuilderTest, build) {
  WingBuilder builder(&containers, vortices);
  builder.AddWing(wing).Build();

  ASSERT_EQ(1, containers.size());
  auto& container = containers[0];
  EXPECT_EQ(2, container.rows());
  EXPECT_EQ(6, container.cols());

  EXPECT_EQ(0, container.at(0, 0).nodes()[0].x());
  EXPECT_EQ(-3, container.at(0, 0).nodes()[0].y());
  EXPECT_EQ(0, container.at(0, 0).nodes()[0].z());
}

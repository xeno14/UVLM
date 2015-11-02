#include <gtest/gtest.h>

#include "../proto/uvlm.pb.h"
#include "proto_adaptor.h"


class Snapshot2ToContainersTest : public ::testing::Test {
  protected:
    virtual void SetUp() {
      snapshot2.Clear();
      for (std::size_t i=0; i<2 ;i++) {
        auto* shape = snapshot2.add_container_shapes();
        shape->set_id(i);
        shape->set_rows(2);
        shape->set_cols(1);
      }
      for (int i=0; i<4; i++) {
        auto* v = snapshot2.add_vortices();
        v->set_gamma(i);
      }
      containers.clear();
      vortices = UVLM::Snapshot2ToContainers(&containers, snapshot2);
    }
    std::vector<UVLM::VortexContainer> containers;
    UVLM::proto::Snapshot2 snapshot2;
    std::shared_ptr<std::vector<UVLM::VortexRing>> vortices;
};

TEST_F(Snapshot2ToContainersTest, vortices) {
  EXPECT_DOUBLE_EQ(0, vortices->at(0).gamma());
  EXPECT_DOUBLE_EQ(1, vortices->at(1).gamma());
  EXPECT_DOUBLE_EQ(2, vortices->at(2).gamma());
  EXPECT_DOUBLE_EQ(3, vortices->at(3).gamma());
}

TEST_F(Snapshot2ToContainersTest, containers) {
  ASSERT_EQ(2, containers.size());
  EXPECT_EQ(2, containers[0].size());
  EXPECT_EQ(2, containers[0].rows());
  EXPECT_EQ(1, containers[0].cols());
  EXPECT_EQ(2, containers[1].size());
  EXPECT_EQ(2, containers[1].rows());
  EXPECT_EQ(1, containers[1].cols());
  
  EXPECT_DOUBLE_EQ(0, containers[0][0].gamma());
  EXPECT_DOUBLE_EQ(1, containers[0][1].gamma());
  EXPECT_DOUBLE_EQ(2, containers[1][0].gamma());
  EXPECT_DOUBLE_EQ(3, containers[1][1].gamma());
}

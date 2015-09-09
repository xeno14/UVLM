#include <gtest/gtest.h>

#include "vortex_container.h"

using UVLM::VortexRing;
using UVLM::VortexContainer;

class VortexContainerTest : public ::testing::Test {
 protected:
  const std::size_t rows, cols;
  VortexContainerTest() : rows(2), cols(3) {}
  virtual void SetUp() {
    // 0 1 2 | 6  7  8 | . . . | . . . 
    // 3 4 5 | 9 10 11 | . . . | . . .
    vortices = std::make_shared<typename decltype(vortices)::element_type>();
    int val = 0;
    for (std::size_t id = 0; id < 4; id++) {
      for (std::size_t i=0; i<rows; i++) {
        for (std::size_t j=0; j<cols; j++) {
          vortices->emplace_back();
          vortices->rbegin()->set_gamma(val++);
        }
      }
    }
  }
  virtual void TearDown() {
    vortices.reset();
  }
  std::shared_ptr<std::vector<VortexRing>> vortices;
  VortexContainer v;
};

TEST_F(VortexContainerTest, constract) {
  VortexContainer v(vortices, rows, cols, 1);
  EXPECT_EQ(2, v.rows());
  EXPECT_EQ(3, v.cols());
  EXPECT_EQ(1, v.id());
  EXPECT_EQ(6, v.Index(0));
  EXPECT_EQ(7, v.Index(0, 1));
}

TEST_F(VortexContainerTest, index) {
  VortexContainer v(vortices, rows, cols, 1);
  EXPECT_EQ(6, v.Index(0));
  EXPECT_EQ(7, v.Index(1));
  EXPECT_EQ(8, v.Index(2));
  EXPECT_EQ(6, v.Index(0, 0));
  EXPECT_EQ(7, v.Index(0, 1));
  EXPECT_EQ(8, v.Index(0, 2));
  EXPECT_EQ(9, v.Index(1, 0));
  EXPECT_EQ(10, v.Index(1, 1));
  EXPECT_EQ(11, v.Index(1, 2));
}

TEST_F(VortexContainerTest, index2) {
  VortexContainer v(vortices, rows, cols, 2);
  EXPECT_EQ(12, v.Index(0));
  EXPECT_EQ(13, v.Index(1));
  EXPECT_EQ(14, v.Index(2));
  EXPECT_EQ(12, v.Index(0, 0));
  EXPECT_EQ(13, v.Index(0, 1));
  EXPECT_EQ(16, v.Index(1, 1));
  EXPECT_EQ(17, v.Index(1, 2));
}

TEST_F(VortexContainerTest, access) {
  VortexContainer v(vortices, rows, cols, 1); 
  EXPECT_DOUBLE_EQ( 6, v[0].gamma());
  EXPECT_DOUBLE_EQ( 7, v[1].gamma());
  EXPECT_DOUBLE_EQ( 8, v[2].gamma());
  EXPECT_DOUBLE_EQ( 9, v[3].gamma());
  EXPECT_DOUBLE_EQ(10, v[4].gamma());
  EXPECT_DOUBLE_EQ(11, v[5].gamma());
  EXPECT_DOUBLE_EQ( 6, v.at(0, 0).gamma());
  EXPECT_DOUBLE_EQ( 7, v.at(0, 1).gamma());
  EXPECT_DOUBLE_EQ( 8, v.at(0, 2).gamma());
  EXPECT_DOUBLE_EQ( 9, v.at(1, 0).gamma());
  EXPECT_DOUBLE_EQ(10, v.at(1, 1).gamma());
  EXPECT_DOUBLE_EQ(11, v.at(1, 2).gamma());
}

TEST_F(VortexContainerTest, access2) {
  VortexContainer v(vortices, rows, cols, 2); 
  EXPECT_EQ(12, v.Offset());
  EXPECT_DOUBLE_EQ(12, v[0].gamma());
  EXPECT_DOUBLE_EQ(13, v[1].gamma());
  EXPECT_DOUBLE_EQ(12, v.at(0, 0).gamma());
  EXPECT_DOUBLE_EQ(13, v.at(0, 1).gamma());
}

TEST_F(VortexContainerTest, shared) {
  VortexContainer v(vortices, rows, cols, 0); 
  EXPECT_DOUBLE_EQ(0, v[0].gamma());
  v[0].set_gamma(10);
  EXPECT_DOUBLE_EQ(10, vortices->at(0).gamma());

  EXPECT_EQ(&(v[0]), &(vortices->at(0)));
  EXPECT_EQ(vortices->begin(), v.begin());

  EXPECT_DOUBLE_EQ(1, v[1].gamma());
  vortices->at(1).set_gamma(10);
  EXPECT_DOUBLE_EQ(10, v[1].gamma());
}

TEST_F(VortexContainerTest, iterator) {
  VortexContainer v(vortices, rows, cols, 1); 
  auto first = v.begin();
  auto last = v.end();
  EXPECT_FALSE(first == last);
  EXPECT_TRUE(first != last);
  EXPECT_DOUBLE_EQ(6, first->gamma()); ++first;
  EXPECT_DOUBLE_EQ(7, first->gamma()); ++first;
  EXPECT_DOUBLE_EQ(8, first->gamma()); ++first;
  EXPECT_DOUBLE_EQ(9, first->gamma()); ++first;
  EXPECT_DOUBLE_EQ(10, first->gamma()); ++first;
  EXPECT_DOUBLE_EQ(11, first->gamma()); ++first;
  EXPECT_TRUE(first == last);
  EXPECT_FALSE(first != last);
}

TEST_F(VortexContainerTest, edge_iterator) {
  VortexContainer v(vortices, rows, cols, 1); 
  auto first = v.edge_begin();
  auto last = v.edge_end();
  EXPECT_FALSE(first == last);
  EXPECT_TRUE(first != last);
  EXPECT_DOUBLE_EQ(9, first->gamma()); ++first;
  EXPECT_DOUBLE_EQ(10, first->gamma()); ++first;
  EXPECT_DOUBLE_EQ(11, first->gamma()); ++first;
  EXPECT_TRUE(first == last);
  EXPECT_FALSE(first != last);
}

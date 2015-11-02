#include <gtest/gtest.h>

#include "vortex_container.h"

using UVLM::VortexRing;
using UVLM::VortexContainer;

class VortexContainerTest : public ::testing::Test {
 public:
   const std::size_t kRows;
   const std::size_t kCols;
   const std::size_t kIdMax;
   const std::size_t kTotalSize;
   const double kChord, kSpan;
  protected:
   VortexContainerTest()
       : kRows(2),
         kCols(3),
         kIdMax(4),
         kTotalSize(kRows * kCols * kIdMax),
         kChord(1),
         kSpan(2) {}
   virtual void SetUp() {
     // 2*3 コンテナが4つ
     // ↓数字はgammaの値
     // 0 1 2 | 6  7  8 | . . . | . . .
     // 3 4 5 | 9 10 11 | . . . | . . .
     vortices = std::make_shared<typename decltype(vortices)::element_type>();
     int val = 0;
     for (std::size_t id = 0; id < kIdMax; id++) {
       for (std::size_t i = 0; i < kRows; i++) {
         for (std::size_t j = 0; j < kCols; j++) {
           vortices->emplace_back();
           vortices->rbegin()->set_gamma(val++);
         }
       }
     }
     for (std::size_t id = 0; id < kIdMax; id++) {
       containers.emplace_back(vortices, kRows, kCols, kIdMax, kChord, kSpan);
     }
   }
   virtual void TearDown() {
     vortices.reset();
     containers.clear();
   }
   std::shared_ptr<std::vector<VortexRing>> vortices;
   VortexContainer v;
   std::vector<VortexContainer> containers;
};

TEST_F(VortexContainerTest, constract) {
  VortexContainer v(vortices, kRows, kCols, 1, 1, 1);
  EXPECT_EQ(2, v.rows());
  EXPECT_EQ(3, v.cols());
  EXPECT_EQ(1, v.id());
  EXPECT_EQ(6, v.Index(0));
  EXPECT_EQ(7, v.Index(0, 1));
}

TEST_F(VortexContainerTest, index) {
  VortexContainer v(vortices, kRows, kCols, 1, 1, 1);
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
  VortexContainer v(vortices, kRows, kCols, 2, 1, 1);
  EXPECT_EQ(12, v.Index(0));
  EXPECT_EQ(13, v.Index(1));
  EXPECT_EQ(14, v.Index(2));
  EXPECT_EQ(12, v.Index(0, 0));
  EXPECT_EQ(13, v.Index(0, 1));
  EXPECT_EQ(16, v.Index(1, 1));
  EXPECT_EQ(17, v.Index(1, 2));
}

TEST_F(VortexContainerTest, access) {
  VortexContainer v(vortices, kRows, kCols, 1, 1, 1);
  EXPECT_DOUBLE_EQ(6, v[0].gamma());
  EXPECT_DOUBLE_EQ(7, v[1].gamma());
  EXPECT_DOUBLE_EQ(8, v[2].gamma());
  EXPECT_DOUBLE_EQ(9, v[3].gamma());
  EXPECT_DOUBLE_EQ(10, v[4].gamma());
  EXPECT_DOUBLE_EQ(11, v[5].gamma());
  EXPECT_DOUBLE_EQ(6, v.at(0, 0).gamma());
  EXPECT_DOUBLE_EQ(7, v.at(0, 1).gamma());
  EXPECT_DOUBLE_EQ(8, v.at(0, 2).gamma());
  EXPECT_DOUBLE_EQ(9, v.at(1, 0).gamma());
  EXPECT_DOUBLE_EQ(10, v.at(1, 1).gamma());
  EXPECT_DOUBLE_EQ(11, v.at(1, 2).gamma());
}

TEST_F(VortexContainerTest, access2) {
  VortexContainer v(vortices, kRows, kCols, 2, 1, 1);
  EXPECT_EQ(12, v.Offset());
  EXPECT_DOUBLE_EQ(12, v[0].gamma());
  EXPECT_DOUBLE_EQ(13, v[1].gamma());
  EXPECT_DOUBLE_EQ(12, v.at(0, 0).gamma());
  EXPECT_DOUBLE_EQ(13, v.at(0, 1).gamma());
}

TEST_F(VortexContainerTest, shared) {
  VortexContainer v(vortices, kRows, kCols, 0, 1, 1);
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
  VortexContainer v(vortices, kRows, kCols, 1, 1, 1);
  auto first = v.begin();
  auto last = v.end();
  EXPECT_FALSE(first == last);
  EXPECT_TRUE(first != last);
  EXPECT_DOUBLE_EQ(6, first->gamma());
  ++first;
  EXPECT_DOUBLE_EQ(7, first->gamma());
  ++first;
  EXPECT_DOUBLE_EQ(8, first->gamma());
  ++first;
  EXPECT_DOUBLE_EQ(9, first->gamma());
  ++first;
  EXPECT_DOUBLE_EQ(10, first->gamma());
  ++first;
  EXPECT_DOUBLE_EQ(11, first->gamma());
  ++first;
  EXPECT_TRUE(first == last);
  EXPECT_FALSE(first != last);
}

TEST_F(VortexContainerTest, edge_iterator) {
  VortexContainer v(vortices, kRows, kCols, 1, 1, 1);
  auto first = v.edge_begin();
  auto last = v.edge_end();
  EXPECT_FALSE(first == last);
  EXPECT_TRUE(first != last);
  EXPECT_DOUBLE_EQ(9, first->gamma());
  ++first;
  EXPECT_DOUBLE_EQ(10, first->gamma());
  ++first;
  EXPECT_DOUBLE_EQ(11, first->gamma());
  ++first;
  EXPECT_TRUE(first == last);
  EXPECT_FALSE(first != last);
}

TEST_F(VortexContainerTest, TotalSize) {
  ASSERT_EQ(kTotalSize, vortices->size());
  EXPECT_EQ(kTotalSize,
            UVLM::CountTotalSize(containers.begin(), containers.end()));

  int n = 100;
  while (n--) vortices->emplace_back();
  ASSERT_NE(kTotalSize, vortices->size());
  EXPECT_EQ(kTotalSize, CountTotalSize(containers.begin(), containers.end()));
}

TEST_F(VortexContainerTest, CopyContainers) {
  std::vector<VortexContainer> copied(containers.size());
  CopyContainers(containers.begin(), containers.end(), copied.begin());

  EXPECT_NE(copied.begin()->vortices().get(), vortices.get());
  for (std::size_t i = 0; i < containers.size(); i++) {
    EXPECT_TRUE(containers[i].ShapeEquals(copied[i]));
    EXPECT_DOUBLE_EQ(containers[i].chord(), copied[i].chord());
    EXPECT_DOUBLE_EQ(containers[i].span(), copied[i].span());
  }
}

TEST_F(VortexContainerTest, Grad) {
  
}

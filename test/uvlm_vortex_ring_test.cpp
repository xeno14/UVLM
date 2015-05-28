#include <gtest/gtest.h>

#include "uvlm_vortex_ring.h"

using UVLM::UVLMVortexRing;

class InitWingTest : public ::testing::Test {
 protected:
  InitWingTest() : dx(0.1) {
    for (std::size_t i=0; i<=rows; i++) {
      for (std::size_t j=0; j<=cols; j++) {
        pos.emplace_back(i*dx, j*dx, 0);
      }
    }
  }

  virtual void SetUp() {}
  virtual void TearDown() {}
  std::vector<Eigen::Vector3d> pos;
  const double dx;
  const std::size_t cols = 10;
  const std::size_t rows = 5;
};

TEST_F(InitWingTest, test1) {
  UVLMVortexRing rings;
  rings.InitWing(pos, cols);

  const auto& v = rings.bound_vortices();
  for (std::size_t i = 0; i < rows; i++) {
    // y>=0
    for (std::size_t j = 0; j < cols; j++) {
      std::size_t idx = j + i * (2 * cols);
      EXPECT_DOUBLE_EQ(i * dx, v[idx].nodes()[0].x());
      EXPECT_DOUBLE_EQ(j * dx, v[idx].nodes()[0].y());
      EXPECT_DOUBLE_EQ(0, v[idx].nodes()[0].z());
      EXPECT_DOUBLE_EQ((i + 1) * dx, v[idx].nodes()[1].x());
      EXPECT_DOUBLE_EQ(j * dx, v[idx].nodes()[1].y());
      EXPECT_DOUBLE_EQ(0, v[idx].nodes()[1].z());
      EXPECT_DOUBLE_EQ((i + 1) * dx, v[idx].nodes()[2].x());
      EXPECT_DOUBLE_EQ((j + 1) * dx, v[idx].nodes()[2].y());
      EXPECT_DOUBLE_EQ(0, v[idx].nodes()[2].z());
      EXPECT_DOUBLE_EQ(i * dx, v[idx].nodes()[3].x());
      EXPECT_DOUBLE_EQ((j + 1) * dx, v[idx].nodes()[3].y());
      EXPECT_DOUBLE_EQ(0, v[idx].nodes()[3].z());
    }
    // y<=0
    for (std::size_t j=0; j<cols; j++) {
      std::size_t idx = j + cols + i*(2*cols);
      EXPECT_DOUBLE_EQ(i * dx, v[idx].nodes()[3].x());
      EXPECT_DOUBLE_EQ(-1.*j * dx, v[idx].nodes()[3].y());
      EXPECT_DOUBLE_EQ(0, v[idx].nodes()[3].z());
      EXPECT_DOUBLE_EQ((i + 1) * dx, v[idx].nodes()[2].x());
      EXPECT_DOUBLE_EQ(-1.*j * dx, v[idx].nodes()[2].y());
      EXPECT_DOUBLE_EQ(0, v[idx].nodes()[1].z());
      EXPECT_DOUBLE_EQ((i + 1) * dx, v[idx].nodes()[1].x());
      EXPECT_DOUBLE_EQ(-1.*(j + 1) * dx, v[idx].nodes()[1].y());
      EXPECT_DOUBLE_EQ(0, v[idx].nodes()[1].z());
      EXPECT_DOUBLE_EQ(i * dx, v[idx].nodes()[0].x());
      EXPECT_DOUBLE_EQ(-1.*(j + 1) * dx, v[idx].nodes()[0].y());
      EXPECT_DOUBLE_EQ(0, v[idx].nodes()[0].z());
    }
  }
}

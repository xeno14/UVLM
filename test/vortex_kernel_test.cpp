#include <gtest/gtest.h>

#include "vortex_kernel.h"
#include "testutil.h"

using UVLM::vortex_kernel::CutOffKernel;

class CutOffTest : public ::testing::Test {};

TEST_F(CutOffTest, cutoff) {
  CutOffKernel kernel(1e-6);
  auto v = kernel.Induce(
      {0, 0, 0},
      {0, 0, -1},
      {0, 0, 1}, 1);
  EXPECT_VECTOR3D_EQ(0, 0, 0, v);
}

TEST_F(CutOffTest, infinite_length) {
  // Sufficiently long line segment induces velocity of G /2 pi r
  CutOffKernel kernel(1e-6);
}

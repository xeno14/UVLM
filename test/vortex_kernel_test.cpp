#include <gtest/gtest.h>

#include "vortex_kernel.h"
#include "testutil.h"

using UVLM::vortex_kernel::CutOffKernel;
using UVLM::vortex_kernel::RosenheadMooreKernel;

const double EPS = 1e-10;

class CutOffTest : public ::testing::Test {
 protected:
  CutOffTest() : kernel(1e-5) {}
  CutOffKernel kernel;
};

TEST_F(CutOffTest, cutoff) {
  auto v = kernel.Induce({0, 0, 0}, {0, 0, -1}, {0, 0, 1}, 1);
  EXPECT_VECTOR3D_EQ(0, 0, 0, v);
}

TEST_F(CutOffTest, infinite_length) {
  // Sufficiently long line segment induces velocity of G /2 pi r

  const double r = 3;
  const double g = 2;
  EXPECT_VECTOR3D_NEAR(0, g / 2 / M_PI / r, 0,
                       kernel.Induce({r, 0, 0}, {0, 0, -1e10}, {0, 0, 1e10}, g),
                       EPS);
}

TEST(RosenheadMooreTest, biotsavart) {
  RosenheadMooreKernel kernel(0);  // same as normal biot-savart

  const double r = 3;
  const double g = 2;
  EXPECT_VECTOR3D_NEAR(0, g / 2 / M_PI / r, 0,
                       kernel.Induce({r, 0, 0}, {0, 0, -1e10}, {0, 0, 1e10}, g),
                       EPS);
}

TEST(RosenheadMooreTest, small_delta) {
  RosenheadMooreKernel kernel(1e-6);

  const double r = 3;
  const double g = 2;
  EXPECT_VECTOR3D_NEAR(0, g / 2 / M_PI / r, 0,
                       kernel.Induce({r, 0, 0}, {0, 0, -1e10}, {0, 0, 1e10}, g),
                       EPS);
}

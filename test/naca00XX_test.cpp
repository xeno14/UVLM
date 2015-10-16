#include <gtest/gtest.h>

#include "testutil.h"
#include "wing/naca00XX.h"

#include <cmath>
#include <fstream>
#include <utility>
#include <vector>
#include <cstdlib>

class NACA0012Test : public ::testing::Test {
 protected:
  NACA0012Test() {
    std::ifstream ifs(UVLM_PROJECT_SOURCE_DIR "/test/wing/data/NACA0012.dat");
    if (!ifs) std::exit(EXIT_FAILURE);
    double x, z;
    ifs >> x >> z;
    data.emplace_back(x, z);
  }
  virtual void SetUp() {}

  virtual void TearDown() {}
  std::vector<std::pair<double, double>> data;
};

TEST_F(NACA0012Test, Defference) {
  double err = 0;
  for (const auto& p : data) {
    double x = p.first;
    double z = p.second;
    err += std::pow(z - UVLM::wing::NACA00XX(x, 1.0, 12), 2);
  }
  // 1e-10 > err
  EXPECT_GT(1e-10, err);
}

TEST_F(NACA0012Test, Generate) {
  ::UVLM::proto::Wing wing;
  ::UVLM::wing::NACA00XXGenerator(12, 2., 2., 2, 2).Generate(&wing);
  const double z1 = UVLM::wing::NACA00XX(1., 2., 12);
  const double z2 = UVLM::wing::NACA00XX(2., 2., 12);
  EXPECT_VECTOR3D_EQ(0, 0, 0, wing.points(0)); 
  EXPECT_VECTOR3D_EQ(0, 1, 0, wing.points(1)); 
  EXPECT_VECTOR3D_EQ(0, 2, 0, wing.points(2)); 
  EXPECT_VECTOR3D_EQ(1, 0, z1, wing.points(3)); 
  EXPECT_VECTOR3D_EQ(1, 1, z1, wing.points(4)); 
  EXPECT_VECTOR3D_EQ(1, 2, z1, wing.points(5)); 
  EXPECT_VECTOR3D_EQ(2, 0, z2, wing.points(6)); 
  EXPECT_VECTOR3D_EQ(2, 1, z2, wing.points(7)); 
  EXPECT_VECTOR3D_EQ(2, 2, z2, wing.points(8)); 
}

TEST_F(NACA0012Test, GenerateRect) {
  ::UVLM::proto::Wing wing;
  ::UVLM::wing::NACA00XXGenerator(0, 2., 2., 2, 2).Generate(&wing);
  EXPECT_VECTOR3D_EQ(0, 0, 0, wing.points(0)); 
  EXPECT_VECTOR3D_EQ(0, 1, 0, wing.points(1)); 
  EXPECT_VECTOR3D_EQ(0, 2, 0, wing.points(2)); 
  EXPECT_VECTOR3D_EQ(1, 0, 0, wing.points(3)); 
  EXPECT_VECTOR3D_EQ(1, 1, 0, wing.points(4)); 
  EXPECT_VECTOR3D_EQ(1, 2, 0, wing.points(5)); 
  EXPECT_VECTOR3D_EQ(2, 0, 0, wing.points(6)); 
  EXPECT_VECTOR3D_EQ(2, 1, 0, wing.points(7)); 
  EXPECT_VECTOR3D_EQ(2, 2, 0, wing.points(8)); 
}

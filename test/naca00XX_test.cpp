#include <gtest/gtest.h>

#include "wing/naca00XX.h"

#include <cmath>
#include <fstream>
#include <utility>
#include <vector>

#include <iostream>
class NACA0012Test : public ::testing::Test {
 protected:
  NACA0012Test() {
    std::ifstream ifs("../../test/wing/data/NACA0012.dat");
    while (!ifs.eof()) {
      double x, z;
      ifs >> x >> z;
      data.emplace_back(x, z);
    }
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

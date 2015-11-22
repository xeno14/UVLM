#include <gtest/gtest.h>

#include "testutil.h"
#include "wing/naca4digit.h"

#include <cmath>
#include <fstream>
#include <utility>
#include <vector>
#include <cstdlib>

class NACA4digitTest : public ::testing::Test {
 protected:
  // Returns mean camber
  auto GetWingData(const char* wingdat) {
    std::vector<std::pair<double, double>> data;
    std::string datapath = UVLM_PROJECT_SOURCE_DIR "/test/wing/data/";
    datapath += wingdat;
    std::ifstream ifs(datapath.c_str());
    if (!ifs) std::exit(EXIT_FAILURE);
    double x, y;
    while(!ifs.eof()) {
      ifs >> x >> y;
      data.emplace_back(x, y);
    }
    return data;
  }
  const double EPS = 1e-6;
};

TEST_F(NACA4digitTest, NACA62XX) {
  auto data = GetWingData("NACA62XX.dat");
  double err = 0;
  for (const auto& p : data) {
    double x = p.first;
    double y = p.second;
    double yc = UVLM::wing::NACA4digit(x, 100, 62);
    std::cerr << x <<" " << y << " " <<yc<<std::endl; 
    err += std::pow(y - yc, 2);
  }
  err /= data.size();
  EXPECT_GT(1e-6, err);
}

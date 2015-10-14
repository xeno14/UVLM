#include <gtest/gtest.h>

#include "testutil.h"
#include "wing/rect.h"

TEST(TransformChebyshev, uni) {
  std::vector<double> x;
  UVLM::wing::TransformChebyshev(&x, 5, 0, 1);
  EXPECT_DOUBLE_EQ(0, x[0]);
  EXPECT_DOUBLE_EQ(cos(-M_PI / 8 * 3), x[1]);
  EXPECT_DOUBLE_EQ(cos(-M_PI / 8 * 2), x[2]);
  EXPECT_DOUBLE_EQ(cos(-M_PI / 8 * 1), x[3]);
  EXPECT_DOUBLE_EQ(1, x[x.size() - 1]);
}

TEST(TransformChebyshev, duo) {
  std::vector<double> x;
  UVLM::wing::TransformChebyshev(&x, 5, 0, 2);
  EXPECT_DOUBLE_EQ(0, x[0]);
  EXPECT_DOUBLE_EQ(2 * cos(-M_PI / 8 * 3), x[1]);
  EXPECT_DOUBLE_EQ(2 * cos(-M_PI / 8 * 2), x[2]);
  EXPECT_DOUBLE_EQ(2 * cos(-M_PI / 8 * 1), x[3]);
  EXPECT_DOUBLE_EQ(2 * 1, x[x.size() - 1]);
}

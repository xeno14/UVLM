#include <gtest/gtest.h>

#include "util.h"
#include <algorithm>

TEST(linspaceTest, plus) {
  auto x = linspace(-1, 1, 3);
  EXPECT_DOUBLE_EQ(-1, x[0]);
  EXPECT_DOUBLE_EQ(0, x[1]);
  EXPECT_DOUBLE_EQ(1, x[2]);
}

TEST(linspaceTest, minus) {
  auto x = linspace(1, -1, 3);
  EXPECT_DOUBLE_EQ(1, x[0]);
  EXPECT_DOUBLE_EQ(0, x[1]);
  EXPECT_DOUBLE_EQ(-1, x[2]);
}

TEST(linspaceTest, cheby) {
  auto theta = linspace(M_PI_2, 0, 10);
  std::vector<double> x(theta.size());
  EXPECT_DOUBLE_EQ(M_PI_2, *(theta.begin()));
  EXPECT_DOUBLE_EQ(0, *(theta.rbegin()));
  std::transform(theta.begin(), theta.end(), x.begin(),
                 [](double t) { return cos(t); });
  x[0] = 0;
  x[1] = 1;
  EXPECT_EQ(0, *(x.begin()));
  EXPECT_EQ(1, *(x.rbegin()));
}

TEST(DoubleLoop, test) {
  auto result = DoubleLoop(2, 3);
  EXPECT_EQ(0, result[0].first);
  EXPECT_EQ(0, result[0].second);
  EXPECT_EQ(0, result[1].first);
  EXPECT_EQ(1, result[1].second);
  EXPECT_EQ(0, result[2].first);
  EXPECT_EQ(2, result[2].second);
  EXPECT_EQ(1, result[3].first);
  EXPECT_EQ(0, result[3].second);
  EXPECT_EQ(1, result[4].first);
  EXPECT_EQ(1, result[4].second);
  EXPECT_EQ(1, result[5].first);
  EXPECT_EQ(2, result[5].second);
}

TEST(JoinTest, test) {
  auto a = {1, 2, 3};
  EXPECT_EQ("1 2 3", UVLM::util::join(" ", a.begin(), a.end()));

  auto b = {"a", "b", "c"};
  EXPECT_EQ("a___b___c", UVLM::util::join("___", b.begin(), b.end()));
}

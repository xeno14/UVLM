#include <gtest/gtest.h>

#include "util.h"

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

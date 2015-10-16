
/**
 * @file testutil.h
 * @brief Add description here
 */
#pragma once

#define EXPECT_VECTOR3D_EQ(X, Y, Z, ACTUAL) \
  EXPECT_DOUBLE_EQ(X, ACTUAL.x()); \
  EXPECT_DOUBLE_EQ(Y, ACTUAL.y()); \
  EXPECT_DOUBLE_EQ(Z, ACTUAL.z())

#define EXPECT_VECTOR3D_NEAR(X, Y, Z, ACTUAL, ERROR) \
  EXPECT_NEAR(X, ACTUAL.x(), ERROR); \
  EXPECT_NEAR(Y, ACTUAL.y(), ERROR); \
  EXPECT_NEAR(Z, ACTUAL.z(), ERROR)

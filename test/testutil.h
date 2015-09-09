
/**
 * @file testutil.h
 * @brief Add description here
 */
#pragma once

#define EXPECT_VECTOR3D_EQ(X, Y, Z, ACTUAL) \
  EXPECT_DOUBLE_EQ(X, ACTUAL.x()); \
  EXPECT_DOUBLE_EQ(Y, ACTUAL.y()); \
  EXPECT_DOUBLE_EQ(Z, ACTUAL.z())

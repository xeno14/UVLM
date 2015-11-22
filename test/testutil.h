
/**
 * @file testutil.h
 * @brief Add description here
 */
#pragma once

#include "vortex.h"

#define EXPECT_VECTOR3D_EQ(X, Y, Z, ACTUAL) \
  EXPECT_DOUBLE_EQ(X, ACTUAL.x()); \
  EXPECT_DOUBLE_EQ(Y, ACTUAL.y()); \
  EXPECT_DOUBLE_EQ(Z, ACTUAL.z())

#define EXPECT_VECTOR3D_NEAR(X, Y, Z, ACTUAL, ERROR) \
  EXPECT_NEAR(X, ACTUAL.x(), ERROR); \
  EXPECT_NEAR(Y, ACTUAL.y(), ERROR); \
  EXPECT_NEAR(Z, ACTUAL.z(), ERROR)

/**
 * @brief Square vortex ring panel shareing a corner with origin
 */
inline UVLM::VortexRing GetSquareRing(double l) {
  UVLM::VortexRing res;
  res.PushNode(Vector3d(0, 0, 0))
     .PushNode(Vector3d(l, 0, 0))
     .PushNode(Vector3d(l, l, 0))
     .PushNode(Vector3d(0, l, 0));
  return res;
}

/**
 * @brief Square vortex ring panel specifing center
 * @param l length of each edge
 * @param x0, y0 center of panel
 */
inline UVLM::VortexRing GetSquareRing(double l, double x0, double y0) {
  UVLM::VortexRing res;
  x0 -= l/2;
  y0 -= l/2;
  res.PushNode(Vector3d(x0, y0, 0))
     .PushNode(Vector3d(x0 + l, y0, 0))
     .PushNode(Vector3d(x0 + l, y0 + l, 0))
     .PushNode(Vector3d(x0, y0 + l, 0));
  return res;
}

#include <gtest/gtest.h>

#include "testutil.h"
#include "calc_impulse.h"

TEST(LambVectorTest, no_ohter_vortex1) {
  std::vector<UVLM::VortexRing> dummy;

  std::vector<Eigen::Vector3d> v_nodes;
  v_nodes.emplace_back(1, 0, 0);
  v_nodes.emplace_back(2, 0, 0);
  v_nodes.emplace_back(3, 0, 0);
  v_nodes.emplace_back(4, 0, 0);

  auto v = GetSquareRing(2);
  v.set_gamma(1);
  auto result = UVLM::calc_impulse::CalcLambVectorOnPanel(
      v, v_nodes, Eigen::Vector3d::Zero(), dummy);
  EXPECT_VECTOR3D_EQ(0, 0, 4, result);
}

TEST(LambVectorTest, no_ohter_vortex2) {
  std::vector<UVLM::VortexRing> dummy;

  std::vector<Eigen::Vector3d> v_nodes;
  v_nodes.emplace_back(1, 0, 0);
  v_nodes.emplace_back(2, 0, 0);
  v_nodes.emplace_back(3, 0, 0);
  v_nodes.emplace_back(4, 0, 0);

  auto v = GetSquareRing(2);
  v.set_gamma(2);
  auto result = UVLM::calc_impulse::CalcLambVectorOnPanel(
      v, v_nodes, Eigen::Vector3d::Zero(), dummy);
  EXPECT_VECTOR3D_EQ(0, 0, 8, result);
}

TEST(LambVectorTest, no_ohter_vortex3) {
  std::vector<UVLM::VortexRing> dummy;

  std::vector<Eigen::Vector3d> v_nodes;
  v_nodes.emplace_back(0, 1, 0);
  v_nodes.emplace_back(0, 2, 0);
  v_nodes.emplace_back(0, 3, 0);
  v_nodes.emplace_back(0, 4, 0);

  auto v = GetSquareRing(2);
  v.set_gamma(1);
  auto result = UVLM::calc_impulse::CalcLambVectorOnPanel(
      v, v_nodes, Eigen::Vector3d::Zero(), dummy);
  EXPECT_VECTOR3D_EQ(0, 0, -4, result);
}

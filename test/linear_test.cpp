#include <gtest/gtest.h>

#include "testutil.h"
#include "linear.h"

#include <cmath>

using UVLM::Morphing;
using UVLM::VortexRing;

namespace {
const double EPS = 1e-10;
}

class CalcMorphingVelocityTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
    m.Clear();
  }

  virtual void TearDown() {

  }
  Morphing m;
};

TEST_F(CalcMorphingVelocityTest, flap) {
  m.set_flap([](double t) { return sin(t);});
  m.set_origin({0, 0, 0});
  std::vector<VortexRing> panels;
  panels.emplace_back(GetSquareRing(1, 0, 1));
  panels.emplace_back(GetSquareRing(1, 0, -2));
  for (auto& v : panels) v.SaveReferenceNode();

  EXPECT_VECTOR3D_EQ(0, 0, 1, panels[0].Normal());
  EXPECT_VECTOR3D_EQ(0, 0, 1, panels[1].Normal());
  
  EXPECT_VECTOR3D_EQ(0, 1, 0, panels[0].ReferenceCentroid());
  EXPECT_VECTOR3D_EQ(0, -2, 0, panels[1].ReferenceCentroid());

  std::vector<double> result(2);
  auto it = UVLM::internal::CalcMorphingVelocityOnWing(
      panels.cbegin(), panels.cend(), m, 0, result.begin());
  EXPECT_EQ(result.end(), it);
  EXPECT_NEAR(1, result[0], EPS);
  EXPECT_NEAR(2, result[1], EPS);

  // TODO test twist
}

class CalcSolveLinearTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
    vortices = std::vector<VortexRing>{GetSquareRing(2, -1, 0),
                                       GetSquareRing(2, 1, 0)};
    for (auto& v : vortices) v.SaveReferenceNode();
  }
  std::vector<VortexRing> vortices;
};

TEST_F(CalcSolveLinearTest, no_wake) {
  // alpha = pi/6
  // no wake
  const double alpha = M_PI/6;
  Morphing m;
  const double U = 2;
  Eigen::Vector3d freestream (U, 0, 0);
  std::vector<VortexRing> vs, dummy;
  VortexRing v;
  const double a = 0.1;
  v.PushNode(0, 0, 0)
      .PushNode(a*cos(alpha), 0, -a*sin(alpha))
      .PushNode(a*cos(alpha), a, -a*sin(alpha))
      .PushNode(0, a, 0);
  v.SaveReferenceNode();
  vs.push_back(v);
  auto result = UVLM::SolveLinearProblem(vs.begin(), vs.end(), dummy.begin(),
                                         dummy.end(), freestream, m, 0);
  double expected = M_SQRT2 * M_PI * a *U * sin(alpha) / 4;
  EXPECT_NEAR(expected, result(0), EPS);
}

TEST_F(CalcSolveLinearTest, no_wake_no_morphing) {
  std::vector<VortexRing> dummy;
  Morphing m;
  auto result = UVLM::SolveLinearProblem(vortices.cbegin(), vortices.cend(),
                                         dummy.cbegin(), dummy.cend(),
                                         Eigen::Vector3d(0, 0, 1), m, 0);
  vortices[0].set_gamma(result(0));
  vortices[1].set_gamma(result(1));
  Eigen::Vector3d v0, v1;
  UVLM::InducedVelocity(&v0, Eigen::Vector3d(-1, 0, 0), vortices.cbegin(),
                        vortices.cend());
  UVLM::InducedVelocity(&v1, Eigen::Vector3d(1, 0, 0), vortices.cbegin(),
                        vortices.cend());
  EXPECT_VECTOR3D_EQ(0, 0, -1, v0);
  EXPECT_VECTOR3D_EQ(0, 0, -1, v1);
}

TEST_F(CalcSolveLinearTest, no_wake_morphing) {
  std::vector<VortexRing> dummy;
  Morphing m;
  m.set_plug([](double t) { return 0.1 * sin(t); });
  auto result = UVLM::SolveLinearProblem(vortices.cbegin(), vortices.cend(),
                                         dummy.cbegin(), dummy.cend(),
                                         Eigen::Vector3d(0, 0, 1), m, 0);
  vortices[0].set_gamma(result(0));
  vortices[1].set_gamma(result(1));
  Eigen::Vector3d v0, v1;
  UVLM::InducedVelocity(&v0, Eigen::Vector3d(-1, 0, 0), vortices.cbegin(),
                        vortices.cend());
  UVLM::InducedVelocity(&v1, Eigen::Vector3d(1, 0, 0), vortices.cbegin(),
                        vortices.cend());
  EXPECT_VECTOR3D_NEAR(0, 0, -1 + 0.1, v0, EPS);
  EXPECT_VECTOR3D_NEAR(0, 0, -1 + 0.1, v1, EPS);
}

TEST_F(CalcSolveLinearTest, wake_morphing) {
  std::vector<VortexRing> wake{GetSquareRing(2, 2, 0)};
  Morphing m;
  m.set_plug([](double t) { return 0.1 * sin(t); });
  auto result = UVLM::SolveLinearProblem(vortices.cbegin(), vortices.cend(),
                                         wake.cbegin(), wake.cend(),
                                         Eigen::Vector3d(0, 0, 1), m, 0);
  vortices[0].set_gamma(result(0));
  vortices[1].set_gamma(result(1));
  Eigen::Vector3d v0, v1;
  UVLM::InducedVelocity(&v0, Eigen::Vector3d(-1, 0, 0), vortices.cbegin(),
                        vortices.cend());
  UVLM::InducedVelocity(&v1, Eigen::Vector3d(1, 0, 0), vortices.cbegin(),
                        vortices.cend());
  Eigen::Vector3d w0, w1;
  UVLM::InducedVelocity(&w0, Eigen::Vector3d(-1, 0, 0), wake.cbegin(),
                        wake.cend());
  UVLM::InducedVelocity(&w1, Eigen::Vector3d(1, 0, 0), wake.cbegin(),
                        wake.cend());
  EXPECT_VECTOR3D_NEAR(0, 0, -1 + 0.1 + w0.z(), v0, EPS);
  EXPECT_VECTOR3D_NEAR(0, 0, -1 + 0.1 + w1.z(), v1, EPS);
}

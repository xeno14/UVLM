#include <gtest/gtest.h>

#include "testutil.h"
#include "linear.h"

#include <cmath>

using UVLM::Morphing;
using UVLM::VortexRing;

namespace {
const double EPS = 1e-6;
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
  EXPECT_NEAR(-1, result[0], EPS);
  EXPECT_NEAR(-2, result[1], EPS);

  // TODO test twist
}

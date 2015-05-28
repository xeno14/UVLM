#include <gtest/gtest.h>

#include "morphing.h"

using UVLM::Morphing;

class MorphingVelocityTest : public ::testing::Test {
 protected:
  const double EPS;
  MorphingVelocityTest() : EPS(1e-6) {
    for (int i=0; i<=10; i++) {
      span.emplace_back(i*0.1, 0, 0);
    }
  }
  virtual void TearDown() {
    morphing.Clear();
  }
  Morphing morphing;
  std::vector<Eigen::Vector3d> span;  // (x,0,0) x in [0,1]で0.1刻み
};

TEST_F(MorphingVelocityTest, plug) {
  morphing.set_plug([](double t) { return sin(t); });
  for (double t=0; t<M_PI; t+=0.1) {
    for (const auto& x0 : span) {
      Eigen::Vector3d v;
      morphing.Velocity(&v, x0, t);
      EXPECT_DOUBLE_EQ(0, v.x());
      EXPECT_DOUBLE_EQ(0, v.y());
      EXPECT_NEAR(cos(t), v.z(), EPS);
    }
  }
}

TEST_F(MorphingVelocityTest, flap) {
  auto phi = [](double t) { return sin(t); };
  morphing.set_flap(phi);
  for (double t = 0; t < M_PI; t += 0.1) {
    for (const auto& x0 : span) {
      Eigen::Vector3d v;
      morphing.Velocity(&v, x0, t);
      EXPECT_DOUBLE_EQ(0, v.x());
      EXPECT_NEAR(-x0.y() * sin(phi(t)) * cos(t), v.y(), EPS);
      EXPECT_NEAR( x0.y() * cos(phi(t)) * cos(t), v.y(), EPS);
    }
  }
}

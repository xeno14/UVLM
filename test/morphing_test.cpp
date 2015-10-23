#include <gtest/gtest.h>

#include "morphing.h"

using UVLM::Morphing;

class MorphingTestBase : public ::testing::Test {
 protected:
  const double EPS;
  const double dt;
  MorphingTestBase() : EPS(1e-6), dt(1e-6) {}
  virtual void TearDown() { morphing.Clear(); }
  Morphing morphing;

  // Morphing functions
  static double Sine(double t) { return sin(t); }
};


class MorphingPerfomeTest : public MorphingTestBase {
 
};

TEST_F(MorphingPerfomeTest, plug) {
  morphing.set_plug(Sine);

  // z(t) = sin(t)
  Eigen::Vector3d x;
  morphing.Perfome(&x, {1, 0, 0}, 0);
  EXPECT_NEAR(1, x.x(), EPS);
  EXPECT_NEAR(0, x.y(), EPS);
  EXPECT_NEAR(0, x.z(), EPS);

  morphing.Perfome(&x, {2, 0, 0}, 0);
  EXPECT_NEAR(2, x.x(), EPS);
  EXPECT_NEAR(0, x.y(), EPS);
  EXPECT_NEAR(0, x.z(), EPS);

  morphing.Perfome(&x, {1, 0, 0}, M_PI/2);
  EXPECT_NEAR(1, x.x(), EPS);
  EXPECT_NEAR(0, x.y(), EPS);
  EXPECT_NEAR(1, x.z(), EPS);
}

TEST_F(MorphingPerfomeTest, plug_origin) {
  morphing.set_plug(Sine);
  Eigen::Vector3d origin {1, 2, 3};
  morphing.set_origin({1, 2, 3});

  // z(t) = sin(t)
  Eigen::Vector3d x;
  morphing.Perfome(&x, {2, 2, 3}, 0);
  EXPECT_NEAR(2, x.x(), EPS);
  EXPECT_NEAR(2, x.y(), EPS);
  EXPECT_NEAR(3, x.z(), EPS);

  morphing.Perfome(&x, {3, 2, 3}, 0);
  EXPECT_NEAR(3, x.x(), EPS);
  EXPECT_NEAR(2, x.y(), EPS);
  EXPECT_NEAR(3, x.z(), EPS);

  morphing.Perfome(&x, {2, 2, 3}, M_PI/2);
  EXPECT_NEAR(2, x.x(), EPS);
  EXPECT_NEAR(2, x.y(), EPS);
  EXPECT_NEAR(4, x.z(), EPS);
}

TEST_F(MorphingPerfomeTest, flap) {
  morphing.set_flap([](double t) { return M_PI_4 * sin(t);});
  morphing.set_origin( {0,0,0});

  // y(t) = r0*sin(phi0 + π/4 * sin(t));
  Eigen::Vector3d x;
  morphing.Perfome(&x, {0, 1, 0}, 0);
  EXPECT_NEAR(0, x.x(), EPS);
  EXPECT_NEAR(1, x.y(), EPS);
  EXPECT_NEAR(0, x.z(), EPS);

  morphing.Perfome(&x, {0, 1, 0}, M_PI_2);
  EXPECT_NEAR(0, x.x(), EPS);
  EXPECT_NEAR(M_SQRT1_2, x.y(), EPS);
  EXPECT_NEAR(-M_SQRT1_2, x.z(), EPS);
}

TEST_F(MorphingPerfomeTest, flap_origin) {
  morphing.set_flap([](double t) { return M_PI_4 * sin(t);});
  morphing.set_origin({0,1,0});

  // y(t) = y0 + r0*sin(phi0 + π/4 * sin(t));
  Eigen::Vector3d x;
  morphing.Perfome(&x, {0, 2, 0}, 0);
  EXPECT_NEAR(0, x.x(), EPS);
  EXPECT_NEAR(2, x.y(), EPS);
  EXPECT_NEAR(0, x.z(), EPS);

  morphing.Perfome(&x, {0, 2, 0}, M_PI_2);
  EXPECT_NEAR(0, x.x(), EPS);
  EXPECT_NEAR(1 + M_SQRT1_2, x.y(), EPS);
  EXPECT_NEAR(-M_SQRT1_2, x.z(), EPS);
}


// MorphingPerfomeTest が通っていれば、plugのみ確かめれば大丈夫なはず
class MorphingVelocityTest : public MorphingTestBase {
 protected:
  const double dx;
  MorphingVelocityTest() : dx(0.1) {
    for (int i=0; i<=10; i++) {
      span.emplace_back(i*dx, 0, 0);
      span_n.emplace_back(-i*dx, 0, 0);
    }
    origin0 << 0, 0, 0;
    origin1 << 0, 1, 0;
  }
  std::vector<Eigen::Vector3d> span;    // (x,0,0) x in [0,1]で0.1刻み
  std::vector<Eigen::Vector3d> span_n;  // span のnegative版
  Eigen::Vector3d origin0, origin1;     // (0,0,0), (0,1,0)
};

TEST_F(MorphingVelocityTest, plug) {
  morphing.set_plug(Sine);
  morphing.set_origin(origin0);
  // z(t)  = sin(t)
  // vz(t) = cos(t)
  Eigen::Vector3d v;
  morphing.Velocity(&v, {1, 0, 0}, 0);
  EXPECT_NEAR(0, v.x(), EPS);
  EXPECT_NEAR(0, v.y(), EPS);
  EXPECT_NEAR(1, v.z(), EPS);

  morphing.Velocity(&v, {2, 0, 0}, 0);
  EXPECT_NEAR(0, v.x(), EPS);
  EXPECT_NEAR(0, v.y(), EPS);
  EXPECT_NEAR(1, v.z(), EPS);

  morphing.Velocity(&v, {1, 0, 0}, M_PI/2);
  EXPECT_NEAR(0, v.x(), EPS);
  EXPECT_NEAR(0, v.y(), EPS);
  EXPECT_NEAR(0, v.z(), EPS);

  morphing.Velocity(&v, {1, 2, 3}, M_PI/4);
  EXPECT_NEAR(0, v.x(), EPS);
  EXPECT_NEAR(0, v.y(), EPS);
  EXPECT_NEAR(M_SQRT1_2, v.z(), EPS);
}

TEST_F(MorphingVelocityTest, plung2) {
  morphing.set_plug(
      [](double t) { return 0.175 * sin(t); });
  Eigen::Vector3d v;
  morphing.Velocity(&v, {100, 2000, 3333}, 0);
  EXPECT_NEAR(0, v.x(), EPS);
  EXPECT_NEAR(0, v.y(), EPS);
  EXPECT_NEAR(0.175, v.z(), EPS);
}

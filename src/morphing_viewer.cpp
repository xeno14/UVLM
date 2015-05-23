/**
 * @file morphing_viewer.cpp
 *
 * @brief 翼の運動の可視化
 *
 * @todo パラメータの設定
 */

#include "morphing.h"

#include <cstdio>
#include <vector>

const double DX = 0.1;

struct SineCurve {
  double A, omega, phi;
  double operator()(double t) { return A * sin(omega * t + phi); }
  SineCurve(double A = 0, double omega = 0, double phi = 0)
      : A(A), omega(omega), phi(phi) {}
};

void InitWing(std::vector<Eigen::Vector3d>* points) {
  points->clear();
  for (int i=0; i<10; i++) {
    for (int j=0; j<20; j++) {
      double x = i * DX;
      double y = j * DX;
      double z = 0;
      points->emplace_back(x, y, z);
    }
  }
}

int main() {
  std::vector<Eigen::Vector3d> wing;
  InitWing(&wing);

  UVLM::Morphing m;
  m.set_flapping(SineCurve(M_PI/6, 1, M_PI/2));
  m.set_pluging(SineCurve(0.5, 1, M_PI/2));

  FILE* fp = stdout;
  for (double t=0; t<20; t+=0.1) {
    for (const auto& x0 : wing) {
      Eigen::Vector3d x;
      m.Perfome(&x, x0, t);
      fprintf(fp, "%e\t%e\t%e\n", x.x(), x.y(), x.z());
      fprintf(fp, "%e\t%e\t%e\n", x.x(), -x.y(), x.z());
    }
    fprintf(fp, "\n\n");
  }

  return 0;
}

/**
 * @file morphing_viewer.cpp
 *
 * @brief 翼の運動の可視化
 *
 * @todo パラメータの設定
 */

#include "morphing.h"

#include <cstdio>
#include <chrono>
#include <thread>
#include <vector>

const double DX = 0.1;

struct SineCurve {
  double A, omega, phi;
  double operator()(double t) { return A * sin(omega * t + phi); }
  SineCurve(double A = 0, double omega = 0, double phi = 0)
      : A(A), omega(omega), phi(phi) {}
};

struct Deformation {
  double s;
  double operator()(const Eigen::Vector3d& x0, double t) {
    return 0.1*sin(t);
  }
};

struct Bend {
  double A, k, omega;
  double operator()(const Eigen::Vector3d& x0, double t) {
    return A * sin(k * x0.y()) * sin(omega * t);
  }
  Bend(double A, double k, double omega) : A(A), k(k), omega(omega) {}
};

void InitWing(std::vector<Eigen::Vector3d>* points) {
  points->clear();
  for (int i=1; i<=10; i++) {
    for (int j=1; j<=40; j++) {
      double x = i * DX;
      double y = -DX/2 + j * DX;
      double z = 0;
      points->emplace_back(x, y, z);
    }
  }
}

int main() {
  std::vector<Eigen::Vector3d> wing;
  InitWing(&wing);

  UVLM::Morphing m;
  m.set_flap(SineCurve(M_PI/6, 1, M_PI/2));
  m.set_plug(SineCurve(0.5, 1, M_PI/2));
  // m.set_twist(Deformation());
  m.set_bend(Bend(0.2, 1, 1));

  FILE* fp = popen("gnuplot", "w");
  fprintf(fp, "set xrange[-1:2]\n");
  fprintf(fp, "set yrange[-3:3]\n");
  fprintf(fp, "set zrange[-2:2]\n");
  for (double t=0; t<1000; t+=0.1) {
    fprintf(fp, "splot \"-\" usi 1:2:3 title \"t=%e\"\n", t);
    // Data
    for (const auto& x0 : wing) {
      Eigen::Vector3d x;
      m.Perfome(&x, x0, t);
      fprintf(fp, "%e\t%e\t%e\n", x.x(), x.y(), x.z());
      fprintf(fp, "%e\t%e\t%e\n", x.x(), -x.y(), x.z());
    }
    fprintf(fp, "end\n");
    fprintf(fp, "\n");
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
  }
  fclose(fp);

  return 0;
}

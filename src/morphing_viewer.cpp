/**
 * @file morphing_viewer.cpp
 *
 * @brief 翼の運動の可視化
 *
 * @todo パラメータの設定
 */

#include "util.h"
#include "morphing.h"
#include "../proto/uvlm.pb.h"

#include <gflags/gflags.h>
#include <cstdio>
#include <chrono>
#include <fstream>
#include <iostream>
#include <thread>
#include <vector>

const double DX = 0.1;


DEFINE_string(wing, "", "Wing proto bin file");


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

void InitWing(std::vector<Eigen::Vector3d>* points, const Eigen::Vector3d& origin) {
  points->clear();

  if (FLAGS_wing.size()) {
    std::ifstream ifs(FLAGS_wing, std::ios::binary);
    CHECK_OPEN(ifs);
    UVLM::proto::Wing wing;
    wing.ParseFromIstream(&ifs);
    for (const auto& point : wing.points()) {
      points->emplace_back(point.x(), point.y(), point.z());
      *points->rbegin() += origin;
      points->emplace_back(point.x(), -point.y(), point.z());
      *points->rbegin() += origin;
    }
  } else {
    for (int i=1; i<=10; i++) {
      for (int j=-20; j<=20; j++) {
        double x = i * DX;
        double y = j * DX;
        double z = 0;
        points->emplace_back(x, y, z);
        *points->rbegin() += origin;
      }
    }
  }
}

int main(int argc, char* argv[]) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  Eigen::Vector3d origin {0, 1, 0};
  const double dt = 0.1;

  std::vector<Eigen::Vector3d> wing;
  InitWing(&wing, origin);

  UVLM::Morphing m;
  m.set_flap(SineCurve(M_PI/6, 1, M_PI/2));
  m.set_plug(SineCurve(0.5, 1, 0));
  // m.set_twist(Deformation());
  m.set_bend(Bend(0.2, 1, 1));

  FILE* fp = popen("gnuplot", "w");
  fprintf(fp, "set xrange[-1:2]\n");
  fprintf(fp, "set yrange[-3:3]\n");
  fprintf(fp, "set zrange[-2:2]\n");
  fprintf(fp, "set xlabel 'x'\n");
  fprintf(fp, "set ylabel 'y'\n");
  fprintf(fp, "set zlabel 'z'\n");
  for (double t=0; t<1000; t+=dt) {
    fprintf(fp, "splot \"-\" usi 1:2:3 title \"t=%e\"\n", t);
    // Data
    for (const auto& x0 : wing) {
      Eigen::Vector3d x;
      m.Perfome(&x, origin, x0, t, dt);
      fprintf(fp, "%e\t%e\t%e\n", x.x(), x.y(), x.z());
    }
    fprintf(fp, "end\n");
    fprintf(fp, "\n");
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
  }
  fclose(fp);

  return 0;
}

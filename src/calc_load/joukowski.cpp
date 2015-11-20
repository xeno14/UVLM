/**
 * @file joukowski.cpp
 * @brief Add description here
 */

#include "joukowski.h"


namespace UVLM {
namespace calc_load {

std::vector<VortexLine> GetLines(const MultipleSheet<Eigen::Vector3d>& pos,
                                 const MultipleSheet<Eigen::Vector3d>& pos_init,
                                 const MultipleSheet<double>& gamma,
                                 std::size_t n) {
  std::vector<VortexLine> res;
  for (std::size_t i = 0; i < gamma.rows(); i++) {
    for (std::size_t j = 0; j < gamma.cols(); j++) {
      const std::vector<Eigen::Vector3d> corner = {
          pos.at(n, i, j), pos.at(n, i, j + 1), pos.at(n, i + 1, j + 1),
          pos.at(n, i + 1, j)};
      const std::vector<Eigen::Vector3d> corner_init = {
          pos_init.at(n, i, j), pos_init.at(n, i, j + 1),
          pos_init.at(n, i + 1, j + 1), pos_init.at(n, i + 1, j)};
      for (std::size_t k = 0; k < corner.size(); k++) {
        if (i == pos.rows() - 1 - 1 && k == 2) continue;  // skip T.E
        res.push_back(VortexLine{
            corner[k], corner[(k + 1) % corner.size()], corner_init[k],
            corner_init[(k + 1) % corner_init.size()], gamma.at(n, i, j)});
      }
    }
  }
  return res;
}

Eigen::Vector3d JoukowskiSteady(
    const std::vector<UVLM::calc_load::VortexLine>& lines,
    const std::vector<Eigen::Vector3d>& U, double t) {
  // steady part
  double Fx = 0, Fy = 0, Fz = 0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : Fx, Fy, Fz)
#endif
  for (std::size_t i = 0; i < lines.size(); i++) {
    const auto& line = lines[i];
    const auto& u = U[i];
    Eigen::Vector3d df = u.cross(line.p1 - line.p0) * line.gamma;
    Fx += df.x();
    Fy += df.y();
    Fz += df.z();
  }
  return Eigen::Vector3d(Fx, Fy, Fz);
}

std::vector<double> CalcPanelArea(const MultipleSheet<Eigen::Vector3d>& pos,
                                  const std::size_t n) {
  std::vector<double> res;
  for (std::size_t i = 0; i < pos.rows() - 1; i++) {
    for (std::size_t j = 0; j < pos.cols() - 1; j++) {
      const double A = ((pos.at(n, i + 1, j) - pos.at(n, i, j))
                            .cross(pos.at(n, i, j + 1) - pos.at(n, i, j)))
                           .norm();
      res.emplace_back(A);
    }
  }
  return res;
}

}  // namespace calc_load
}  // namespace UVLM

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

}  // namespace calc_load
}  // namespace UVLM

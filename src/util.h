/**
 * @file util.h
 * @brief Add description here
 */
#pragma once

#include <cmath>
#include <utility>
#include <vector>

#ifndef CHECK_OPEN
#define CHECK_OPEN(fp)                                                         \
  if (!fp) {                                                                   \
    std::cerr << "Open Error at " << __FILE__ << ":" << __LINE__ << std::endl; \
    std::exit(EXIT_FAILURE);                                                   \
  }
#endif

inline double DegToRad(double deg) {
  return deg / 180.0 * M_PI;
}

inline std::vector<double> linspace(double a, double b, std::size_t N) {
  double dx = (b - a) / (N - 1);
  std::vector<double> res;
  for (std::size_t i = 0; i < N; i++) {
    res.push_back(a + i * dx);
  }
  return res;
}

inline std::vector<std::pair<std::size_t, std::size_t>> DoubleLoop(
    std::size_t imax, std::size_t jmax) {
  std::vector<std::pair<std::size_t, std::size_t>> res;
  for (std::size_t i = 0; i < imax; ++i) {
    for (std::size_t j = 0; j < jmax; ++j) {
      res.emplace_back(i, j);
    }
  }
  return res;
}

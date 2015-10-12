/**
 * @file util.h
 * @brief Add description here
 */
#pragma once

#include <cmath>
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

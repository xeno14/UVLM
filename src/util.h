/**
 * @file util.h
 * @brief Add description here
 */
#pragma once

#include <cmath>

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

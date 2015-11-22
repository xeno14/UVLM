/**
 * @file util.h
 * @brief Add description here
 */
#pragma once


#include <boost/range/iterator_range.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <cmath>
#include <iterator>
#include <type_traits>
#include <utility>
#include <vector>
#include <sstream>

#ifndef CHECK_OPEN
#define CHECK_OPEN(fp)                                                         \
  if (!fp) {                                                                   \
    std::cerr << "Open Error at " << __FILE__ << ":" << __LINE__ << std::endl; \
    std::exit(EXIT_FAILURE);                                                   \
  }
#endif

#define ITERATOR_VALUETYPE_CHECK(type, iterator)                             \
  static_assert(                                                             \
      std::is_same<                                                          \
          type, typename std::iterator_traits<iterator>::value_type>::value, \
      "value_type missmatch: give " #type)

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

namespace UVLM {
namespace util {

template <class InputIterator>
std::string join(const std::string& sep, InputIterator first,
                 InputIterator last) {
  std::stringstream ss;
  while (first != last) {
    ss << *first << sep;
    ++first;
  }
  std::string res(ss.str());
  res.resize(res.size() - sep.size());
  return res;
}

template <class... T>
auto zip(const T&... containers) {
  auto zip_begin =
      boost::make_zip_iterator(boost::make_tuple(std::begin(containers)...));
  auto zip_end =
      boost::make_zip_iterator(boost::make_tuple(std::end(containers)...));
  return boost::make_iterator_range(zip_begin, zip_end);
}

}  // namespace util
}  // namespace UVLM

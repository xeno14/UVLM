/**
 * @file naca00XX_main.cpp
 * @brief Add description here
 */

#include "naca00XX.h"

#include <gflags/gflags.h>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

DEFINE_int32(digit, 12, "Digit of NACA00XX");
DEFINE_double(chord, 1.0, "Chord length (x)");
DEFINE_double(span, 2.0, "Span length (y)");
DEFINE_int32(rows, 10, "Number of rows (x).");
DEFINE_int32(cols, 20, "Number of columns (y).");

std::string GenFilename(const std::string& base_path, int digit,
                        std::size_t rows, std::size_t cols) {
  std::stringstream ss;
  ss << base_path << "/"
     << "naca00" << digit << "_" << cols << "x" << rows;
  return ss.str();
}

int main(int argc, char* argv[]) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  UVLM::wing::NACA00XXGenerator generator(
      FLAGS_digit, FLAGS_chord, FLAGS_span, FLAGS_rows, FLAGS_cols);
  generator.set_verbose(true);

  UVLM::proto::Wing wing;
  generator(&wing);

  const std::string filename(
      GenFilename(".", FLAGS_digit, FLAGS_rows, FLAGS_cols));
  std::ofstream ofs(filename, std::ios::binary);
  if (!ofs) {
    std::cerr << "Unable to open " << filename;
    std::exit(EXIT_FAILURE);
  }
  wing.SerializeToOstream(&ofs);

  return 0;
}
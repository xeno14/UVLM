/**
 * @file naca00XX_main.cpp
 * @brief Add description here
 */

#include "naca00XX.h"
#include "../util.h"

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
DEFINE_string(output, "", "Output file.");

std::string GenFilename(const std::string& base_path, int digit,
                        std::size_t rows, std::size_t cols) {
  if (FLAGS_output.size()) {
    return FLAGS_output;
  }
  std::stringstream ss;
  ss << base_path << "/"
     << "naca00" << digit << "_" << cols << "x" << rows;
  return ss.str();
}

void ToStdOut(const UVLM::proto::Wing& wing) {
  for (auto it = wing.points().begin(); it != wing.points().end(); ++it) {
    std::cout << it->x() << " " << it->y() << " " << it->z() << std::endl;
  }
}

int main(int argc, char* argv[]) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  UVLM::wing::NACA00XXGenerator generator(FLAGS_digit);
      
  generator.set_verbose(true);

  UVLM::proto::Wing wing;
  generator(&wing, FLAGS_chord, FLAGS_span, FLAGS_rows, FLAGS_cols);

  const std::string filename(
      GenFilename(".", FLAGS_digit, FLAGS_rows, FLAGS_cols));
  std::ofstream ofs(filename, std::ios::binary);
  CHECK_OPEN(ofs);

  wing.SerializeToOstream(&ofs);

  return 0;
}

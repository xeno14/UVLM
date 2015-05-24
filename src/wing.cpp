/**
 * @file wing.cpp
 * @brief Add description here
 */

#include "wing.h"

#include <iostream>
#include <gflags/gflags.h>


DEFINE_string(output, "stdout", "Where to output");
DEFINE_double(chord, 1.0, "Chord length");
DEFINE_double(aspect_ratio, 4.0, "Aspect ratio");


void GenNACA0012(std::ostream& os, double chord, double ar) {
  int row = 
}

int main(int argc, char *argv[]) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  std::cout << FLAGS_chord << std::endl;

  return 0;
}

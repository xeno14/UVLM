/**
 * @file main.cpp
 * @brief Add description here
 */

#include "wing.h"

#include <gflags/gflags.h>

#include <iostream>
#include <fstream>
#include <functional>
#include <memory>
#include <string>
#include <vector>

DEFINE_string(output, "stdout", "\"stdout\", \"stderr\" or a file path.");
DEFINE_string(format, "tsv", "tsv, proto");


void OutputTsv(std::ostream* os, const UVLM::proto::Wing& wing);

int main(int argc, char* argv[]) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  if (FLAGS_output.empty()) {
    std::cerr << "Empty output path\n";
    return EXIT_FAILURE;
  }

  UVLM::proto::Wing wing;
  UVLM::wing::Perfome(&wing);

  // std::vector<std::ostream*> gc;
  //
  // std::ostream* os;
  // if (FLAGS_output == "stdout") {
  //   os = &std::cout;
  // } else if (FLAGS_output == "stderr") {
  //   os = &std::cerr;
  // } else {
  //   os = new std::ofstream(FLAGS_output, std::ios::binary);
  //   gc.push_back(os);
  //   if (!(*os)) {
  //     std::cerr << "File open error for " << FLAGS_output << std::endl;
  //     return EXIT_FAILURE;
  //   }
  // }
  //
  //
  // for (auto* p : gc) {
  //   delete p;
  // }
  return 0;
}

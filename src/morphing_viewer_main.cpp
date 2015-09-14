/**
 * @file morphing_viewer_main.cpp
 * @brief Add description here
 */

#include "morphing_viewer.h"


int main(int argc, char* argv[]) {
  using namespace morphing_viewer;
  set_plug([](double t) { return sin(t); });
  set_flap([](double t) { return M_PI/6 * sin(t + M_PI/2); });
  return morphing_viewer::main(argc, argv);
}

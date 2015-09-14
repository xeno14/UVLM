
/**
 * @file morphing_viewer.h
 * @brief Add description here
 */
#pragma once

#include "morphing.h"
#include "../proto/uvlm.pb.h"

#include <vector>

namespace morphing_viewer {

int main(int argc, char* argv[]);

void set_plug(std::function<double(double)> f);
void set_flap(std::function<double(double)> f);
void set_twist(std::function<double(const Eigen::Vector3d&, double)> f);
void set_bend(std::function<double(const Eigen::Vector3d&, double)> f);

}  // namespace morphing_viewer

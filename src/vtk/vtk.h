
/**
 * @file vtk.h
 * @brief Add description here
 */
#pragma once

#include <cstdio>

#include "../../proto/uvlm.pb.h"

namespace UVLM {
namespace vtk {

void WriteSnapshot2(FILE* fp, const proto::Snapshot2& snapshot); 
 
}  // namespace vtk
}  // namespace UVLM

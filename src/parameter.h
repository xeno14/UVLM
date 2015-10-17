
/**
 * @file parameter.h
 * @brief パラメーターの読み取り関係
 */
#pragma once

#include <yaml-cpp/yaml.h>

#ifndef DEFINE_PARAM
#define DEFINE_PARAM(type, name, node)            \
  const type PARAM_##name = node[#name].as<type>()
#endif

#ifndef DEFINE_PARAM_VERBOSE
#define DEFINE_PARAM_VERBOSE(type, name, node)               \
  LOG(INFO) << "set value of " << #name;                     \
  DEFINE_PARAM(type, name, node);                            \
  LOG(INFO) << "@param " << #node << "[\"" << #name << "\"]" \
            << ": " << PARAM_##name
#endif

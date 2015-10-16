
/**
 * @file parameter.h
 * @brief パラメーターの読み取り関係
 */
#pragma once

#include <yaml-cpp/yaml.h>

#ifndef DEFINE_PARAM
#define DEFINE_PARAM(type, name, node)                   \
  if (node[#name].Type() == YAML::NodeType::Null) {      \
    LOG(FATAL) << "Parameter '" << #name << "' is null"; \
  }                                                      \
  const type PARAM_##name = node[#name].as<type>()
#endif

#ifndef DEFINE_PARAM_VERBOSE
#define DEFINE_PARAM_VERBOSE(type, name, node)               \
  DEFINE_PARAM(type, name, node);                            \
  LOG(INFO) << "@param " << #node << "[\"" << #name << "\"]" \
            << ": " << PARAM_##name
#endif

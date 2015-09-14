
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
  const type PARAM_##name = node[#name].as<type>();
#endif

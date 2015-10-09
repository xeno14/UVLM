
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

#include <map>
#include <string>

namespace UVLM {
namespace parameter {
namespace internal {

class Converter {
 public:
  virtual ~Converter() {}
  virtual void Convert(void* ptr, const YAML::Node& node) = 0;
};

class IntConverter : public Converter {
 public:
  virtual void Convert(void* ptr, const YAML::Node& node) {
    *((int*)ptr) = node.as<int>();
  }
};

struct Param {
  void* ptr;
  Converter* converter;
};

class ParamManager {
 private:
  std::map<std::string, Param> params_;

  ~ParamManager() {
    for (auto it = params_.begin(); it != params_.end(); ++it) {
      delete it->second.converter;
    }
  }

 public:
  static ParamManager& instance() {
    static ParamManager p;
    return p;
  }

  void Register(std::string name, void* ptr, Converter* converter) {
    params_[name] = Param{ptr, converter};
  }

  void InitParam(const YAML::Node& node) {
    for (auto it=params_.begin(); it!=params_.end(); ++it) {
      std::string key = it->first;
      auto* value_ptr = it->second.ptr;
      auto* converter = it->second.converter;
      converter->Convert(value_ptr, node);
    }
  }
};

class Registerer {
 public:
  template <class T>
  Registerer(const char* name, T* ptr, Converter* converter) {
    ParamManager::instance().Register(name, (void*)ptr, converter);
  }
};

}  // namespace internal

void InitParam(const YAML::Node& node);

}  // namespace parameter
}  // namespace UVLM

#define DEFINE_PARAM_int(name, default_value)                              \
  int PARAM_##name = default_value;                                        \
  UVLM::parameter::internal::Registerer _UVLM_parameter_Registerer_##name( \
      #name, &PARAM_##name, new UVLM::parameter::internal::IntConverter());

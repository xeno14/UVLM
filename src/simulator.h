
/**
 * @file simulator.h
 * @brief Add description here
 */
#pragma once

#include "../proto/uvlm.pb.h"
#include "vortex.h"
#include "morphing.h"

namespace UVLM {
namespace simulator {
namespace internal {

void CheckReady();

void CreateContainers();
auto ShedProcess();

void MorphingProcess(const double t, const double dt);

void AdvectProcess(const double dt);

void AppendShedProcess(std::vector<std::vector<::UVLM::VortexRing>>* shed);

/**
 * @biref 線形問題を解く
 *
 * @todo morphingが１つしか使えないので複数対応する
 */
void SolveLinearProblem();

}  // namespace internal

void InitSimulator();

void AddWing(const ::UVLM::proto::Wing& wing, const UVLM::Morphing& m);

/**
 * @brief 上流の流速を設定する 
 * @todo time dependency
 */
void SetInlet(double x, double y, double z);

void SetOutputPath(const std::string& path);

void Start(std::size_t steps, const double dt);

}  // simulator
}  // namespace UVLM

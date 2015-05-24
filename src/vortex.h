
/**
 * @file vortex.h
 * @brief Add description here
 */
#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>


namespace UVLM {

double DistanceLineAndPoint(const Eigen::Vector3d& start,
                            const Eigen::Vector3d& end,
                            const Eigen::Vector3d& pos);

void BiotSavartLaw(Eigen::Vector3d* result,
                   const Eigen::Vector3d& start, const Eigen::Vector3d& end,
                   const Eigen::Vector3d& pos);

class VortexFilament {
 public:
  VortexFilament() {}
  VortexFilament(double gamma, const Eigen::Vector3d& start,
                 const Eigen::Vector3d& end)
      : gamma_(gamma), start_(start), end_(end) {}
  VortexFilament(const VortexFilament& v)
      : gamma_(v.gamma_), start_(v.start_), end_(v.end_) {}

  void BiotSavartLaw(Eigen::Vector3d* result, const Eigen::Vector3d& pos) const;

 private:
  double gamma_;
  Eigen::Vector3d start_, end_;
};

}  // namespace UVLM

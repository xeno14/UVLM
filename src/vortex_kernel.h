/**
 * @file vortex_kernel.h
 * @brief Add description here
 */
#pragma once

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>

namespace UVLM {
namespace vortex_kernel {

class VortexKernel {
 public:
  virtual ~VortexKernel() = default;
  virtual Eigen::Vector3d Induce(const Eigen::Vector3d& x,
                                 const Eigen::Vector3d& x1,
                                 const Eigen::Vector3d& x2,
                                 const double gamma) const = 0;
};

/**
 * @brief biot-savart law using cutoff
 *
 * Corresponds to subroutine VORTEX in Katz and Plotkin (p.584, second edition).
 * Mathmatical formulation is written at p.255
 */
class CutOffKernel : public VortexKernel {
 public:
  CutOffKernel(const double l = 1e-10) : cutoff_length_(l) {}

  Eigen::Vector3d Induce(const Eigen::Vector3d& x, const Eigen::Vector3d& x1,
                         const Eigen::Vector3d& x2,
                         const double gamma) const override;

 private:
  const double cutoff_length_;
};

}  // namespace vortex_kernel
}  // namespace UVLM

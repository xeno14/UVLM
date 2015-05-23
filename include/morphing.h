/**
 * @file morphing.h
 * @brief Add description here
 */
#pragma once

// std
#include <cmath>
#include <functional>

// Eigen
#include <Eigen/Core>
#include <Eigen/Geometry>

namespace UVLM {
namespace internal {
inline double DefaultFunc(double) { return 0; }
inline double DefaultFunc2(const Eigen::Vector3d&, double) { return 0; }
}

/**
 * @brief 翼の変形を司るクラス
 */
class Morphing {
 public:
  Morphing();
  ~Morphing() = default;

  /** @brief 変形を行う
   *  @param x  変形後の位置
   *  @param x0 reference frameでの座標
   *  @param t  時刻
   */
  void Perfome(Eigen::Vector3d* x, const Eigen::Vector3d& x0,
               const double t) const;

  /** @brief 回転行列を用意する
   *  
   *  Ghommem (2012) Eq. (15)
   */
  void PrepareMatrix(Eigen::Matrix3d* m, const Eigen::Vector3d& x0,
                     double t) const;

  void set_alpha(double alpha) { alpha_ = alpha; }
  void set_plug(std::function<double(double)> f) { plug_ = f; }
  void set_flap(std::function<double(double)> f) { flap_ = f; }
  void set_twist(std::function<double(const Eigen::Vector3d&, double)> f) {
    twist_ = f;
  }

 private:
  double alpha_;
  std::function<double(double)> plug_;
  std::function<double(double)> flap_;
  std::function<double(const Eigen::Vector3d&, double)> twist_;
};

}  // namespace UVLM

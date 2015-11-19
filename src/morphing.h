/**
 * @file morphing.h
 * @brief Add description here
 */
#pragma once

// std
#include <cmath>
#include <functional>

// Eigen
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>

namespace UVLM {
namespace internal {
inline double DefaultFunc(double) { return 0; }
inline double DefaultFunc2(const Eigen::Vector3d&, double) { return 0; }
// inline void DefaultFunc3(Eigen::Vector3d*, double) {}
}

/**
 * @brief 翼の変形を司るクラス
 */
class Morphing {
 public:
  Morphing();
  Morphing(const Morphing&);
  ~Morphing() = default;

  /** @brief 変形を行う
   *  @param x[out]  変形後の位置
   *  @param x0[in] reference frameでの座標
   *  @param t[in]  時刻
   */
  void Perfome(Eigen::Vector3d* x, const Eigen::Vector3d& x0,
               const double t) const;

  /** @brief 変形速度
   *
   *  変形を行い、数値微分（中心差分）をして速度とする。
   *  @param v[out] 変形速度
   *  @param x0[in] reference frameでの位置
   *  @param t[in] 時刻
   */
  void Velocity(Eigen::Vector3d* v,
                const Eigen::Vector3d& x0, const double t,
                const double dt = 1e-6) const;

  Eigen::Vector3d Velocity(const Eigen::Vector3d& x0, const double t,
                           const double dt = 1e-6) const {
    Eigen::Vector3d res;
    Velocity(&res, x0, t, dt);
    return res;
  }

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
  void set_bend(std::function<double(const Eigen::Vector3d&, double)> f) {
    bend_ = f;
  }

  /** @brief Set position of the root of wing */
  void set_origin(const Eigen::Vector3d& origin) { origin_ = origin; }

  /**
   * @brief Clear functions that do nothing
   */
  void Clear();

 private:
  double alpha_;
  std::function<double(double)> plug_;
  std::function<double(double)> flap_;
  std::function<double(const Eigen::Vector3d&, double)> twist_;
  std::function<double(const Eigen::Vector3d&, double)> bend_;
  // std::function<void(Eigen::Vector3d*, double)> thrust_;

  Eigen::Vector3d origin_;
};

}  // namespace UVLM

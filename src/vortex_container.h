/**
 * @file vortex_container.h
 * @brief Add description here
 */
#pragma once

#include "vortex.h"

#include <memory>
#include <vector>

#include <iostream>

namespace UVLM {

/**
 * Bound Vorticesを保持するコンテナで１つの翼に対応する
 *
 * 渦の本体自体はshared_ptrにvectorで保持し、VortexContainerはその中でindexが
 * 何番目から何番目までかを記録する。
 */
class VortexContainer {
  typedef std::vector<VortexRing> container_t;
  typedef std::shared_ptr<container_t> vortices_ptr_t;

 public:
  VortexContainer() : vortices_(nullptr) {}
  VortexContainer(const vortices_ptr_t& vortices, std::size_t rows,
                  std::size_t cols, std::size_t id, double chord, double span)
      : vortices_(vortices),
        rows_(rows),
        cols_(cols),
        id_(id),
        chord_(chord),
        span_(span) {}
  VortexContainer(const VortexContainer& container)
      : vortices_(container.vortices_),
        rows_(container.rows_),
        cols_(container.cols_),
        id_(container.id_),
        chord_(container.chord_),
        span_(container.span_) {}

  VortexContainer& operator=(const VortexContainer& c) {
    vortices_ = c.vortices_;
    rows_ = c.rows_;
    cols_ = c.cols_;
    id_ = c.id_;
    chord_ = c.chord_;
    span_ = c.span_;
    return *this;
  }

  // vorticesの始点からのオフセット
  std::size_t Offset() const { return id_ * rows_ * cols_; }

  // vorticesのindexへの変換
  std::size_t Index(std::size_t i) const { return Offset() + i; }

  // vorticesのindexへの変換
  std::size_t Index(std::size_t i, std::size_t j) const {
    return Offset() + j + i * cols_;
  }

  VortexRing& operator[](std::size_t i) { return (*vortices_)[Index(i)]; }
  const VortexRing& operator[](std::size_t i) const {
    return (*vortices_)[Index(i)];
  }

  VortexRing& at(std::size_t i, std::size_t j) {
    return (*vortices_)[Index(i, j)];
  }
  const VortexRing& at(std::size_t i, std::size_t j) const {
    return (*vortices_)[Index(i, j)];
  }

  bool ShapeEquals(const VortexContainer& c) const {
    return rows_ == c.rows() && cols_ == c.cols() && id_ == c.id();
  }

  void alloc(std::size_t sz) {
    if (vortices_->size() < sz) {
      vortices_->resize(sz);
    }
  }

  double DeltaGamma(std::size_t i, std::size_t j) const {
    const double g = this->at(i, j).gamma();
    const double g_i = i == 0 ? 0 : this->at(i - 1, j).gamma(); // leading edge
    return g - g_i;
  }

  Eigen::Vector3d Grad(std::size_t i, std::size_t j) const {
    const auto& v = this->at(i, j);
    return v.Tangent() * DeltaGamma(i,j) / v.CalcC();
  }

  const vortices_ptr_t& vortices() const { return vortices_; }
  std::size_t cols() const { return cols_; }
  std::size_t rows() const { return rows_; }
  std::size_t id() const { return id_; }
  double chord() const { return chord_; }
  double span() const { return span_; }

  void set_vortices(const vortices_ptr_t& vortices, std::size_t r,
                    std::size_t c, std::size_t i, double chord, double span) {
    vortices_ = vortices;
    rows_ = r;
    cols_ = c;
    id_ = i;
    chord_ = chord;
    span_ = span;
  }

  std::size_t size() const { return cols_ * rows_; }

  auto begin() { return vortices_->begin() + Index(0); }
  auto begin() const { return vortices_->begin() + Index(0); }
  auto end() { return vortices_->begin() + Index(rows_ * cols_); }
  auto end() const { return vortices_->begin() + Index(rows_ * cols_); }
  auto edge_begin() { return vortices_->begin() + Index(rows_ - 1, 0); }
  auto edge_end() { return end(); }
  auto edge_begin() const { return vortices_->begin() + Index(rows_ - 1, 0); }
  auto edge_end() const { return end(); }
  auto cbegin() const { return vortices_->cbegin() + Index(0); }
  auto cend() const { return vortices_->cbegin() + Index(rows_ * cols_); }

  std::vector<Eigen::Vector3d> DumpPos() const;
  void LoadPos(const std::vector<Eigen::Vector3d>& pos);

 private:
  vortices_ptr_t vortices_;
  std::size_t rows_, cols_, id_;
  double chord_, span_;
};

/**
 * コンテナの集合から要素の総数を計算する
 */
template <class InputIterator>
std::size_t CountTotalSize(InputIterator first, InputIterator last) {
  std::size_t res = 0;
  while (first != last) {
    res += first->size();
    ++first;
  }
  return res;
}

/**
 * コンテナをコピーする
 *
 * 共有された渦も新しく領域を確保してコピーする
 */
template <class InputIterator, class OutputIterator>
void CopyContainers(InputIterator first, InputIterator last,
                    OutputIterator result) {
  auto vortices = first->vortices();
  const auto total_size = CountTotalSize(first, last);

  // 渦をコピーする
  auto copied = std::make_shared<std::vector<UVLM::VortexRing>>(
      vortices->begin(), vortices->begin() + total_size);

  while (first != last) {
    *result = VortexContainer(copied, first->rows(), first->cols(), first->id(),
                             first->chord(), first->span());
    ++first;
    ++result;
  }
}

/**
 * @brief wakeのイテレータを入手する
 * @tparam Range VortexContainerのレンジ
 * @return pair {wake_first, wake_last}
 */
template <class Range>
auto GetWake(const Range& containers) {
  const auto wake_offset =
      CountTotalSize(std::begin(containers), std::end(containers));
  // TODO what if vortices is nullptr?
  auto vortices = std::begin(containers)->vortices();
  return std::make_pair(vortices->begin() + wake_offset, vortices->end());
}

}  // namespace UVLM

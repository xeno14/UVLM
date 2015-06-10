/**
 * @file vortex_container.h
 * @brief Add description here
 */
#pragma once

#include "vortex.h"

#include <memory>
#include <vector>


namespace UVLM {

class VortexContainer {
  typedef std::vector<VortexRing> container_t;
  typedef std::shared_ptr<container_t> vortices_ptr_t;

 public:
  VortexContainer() : vortices_(nullptr) {}
  VortexContainer(const vortices_ptr_t& vortices, std::size_t rows,
                  std::size_t cols, std::size_t id)
      : vortices_(vortices), rows_(rows), cols_(cols), id_(id) {}
  VortexContainer(const VortexContainer& container)
      : vortices_(container.vortices_),
        rows_(container.rows_),
        cols_(container.cols_),
        id_(container.id_) {}

  VortexContainer& operator=(const VortexContainer& container) {
    this->set_vortices(container.vortices(), container.rows(), container.cols(),
                       container.id());
    return *this;
  }

  // vorticesの始点からのオフセット
  std::size_t Offset() const {
    return rows_ * cols_;
  }

  // vorticesのindexへの変換
  std::size_t Index(std::size_t i) const {
    return Offset() + i;
  }

  // vorticesのindexへの変換
  std::size_t Index(std::size_t i, std::size_t j) const {
    return Offset() + j + i * cols_;
  }


  VortexRing& operator[] (std::size_t i) { return (*vortices_)[Index(i)]; }
  const VortexRing& operator[] (std::size_t i) const {
    return (*vortices_)[Index(i)];
  }

  VortexRing& at(std::size_t i, std::size_t j) {
    return (*vortices_)[Index(i, j)];
  }
  const VortexRing& at(std::size_t i, std::size_t j) const {
    return (*vortices_)[Index(i, j)];
  }

  const vortices_ptr_t& vortices() const { return vortices_; }
  std::size_t cols() const { return cols_; }
  std::size_t rows() const { return rows_; }
  std::size_t id() const { return id_; }

  void set_vortices(const vortices_ptr_t& vortices, std::size_t r,
                    std::size_t c, std::size_t i) {
    vortices_ = vortices;
    rows_ = r;
    cols_ = c;
    id_ = i;
  }

  auto begin() {
    return vortices_->begin() + Index(0);
  }
  auto end() {
    return vortices_->begin() + Index(rows_ * cols_);
  }
  auto edge_begin() {
    return vortices_->begin() + Index(rows_ - 1, 0);
  }
  auto edge_end() {
    return end();
  }
  auto cbegin() const {
    return vortices_->cbegin() + Index(0);
  }
  auto cend() const {
    return vortices_->cbegin() + Index(rows_ * cols_);
  }

 private:
  vortices_ptr_t vortices_;
  std::size_t rows_, cols_, id_;

};

}  // UVLM

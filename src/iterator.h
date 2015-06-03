
/**
 * @file iterator.h
 * @brief Add description here
 */
#pragma once

#include <iterator>

namespace UVLM {

/** @brief 配列の一部分を取り出すためのイテレータ
 *  @tparam RandomAccessList ランダムアクセス可能なリスト
 */
template <class RandomAccessList>
class Iterator : public std::iterator<std::forward_iterator_tag,
                                      typename RandomAccessList::value_type> {
 public:
  Iterator(RandomAccessList& list, std::size_t idx, std::size_t diff) : list_(list), idx_(idx), diff_(diff) {}
  Iterator(const Iterator& it) : list_(it.list_), idx_(it.idx_), diff_(it.diff_) {}
  Iterator& operator++() {
    idx_ += diff_;
    return *this;
  }
  Iterator operator++(int) {
    Iterator it(*this);
    operator++();
    return it;
  }
  bool operator==(const Iterator& it) const { return idx_ == it.idx_; }
  bool operator!=(const Iterator& it) const { return idx_ != it.idx_; }
  auto& operator*() { return list_[idx_]; }

 private:
  RandomAccessList list_;
  std::size_t idx_, diff_;
};

}  // namespace UVLM

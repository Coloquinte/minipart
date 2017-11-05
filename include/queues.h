// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

/*
 * This file implements priority queues with key updates, used in local search algorithms
 * Most use the same few tricks:
 *    - bucket-based: elements with the same priorities are put in an array
 *    - lazy: elements are pushed each time they are updated; consistency with their priority is checked at pop
 *    - LIFO tie resolution: priority ties are resolved in last-in first-out order
 */

#pragma once

#include "inc_bipart.h"

namespace minipart {

template <bool FIFO>
class BasicQueue;

template <>
class BasicQueue<false> {
 // Basically a stack
 public:
  BasicQueue() {}
  BasicQueue(IncBipart &) {}

  bool empty() const { return d_.empty(); }
  Node top() const { return d_.back(); }
  Node pop() {
    Node n = d_.back();
    d_.pop_back();
    return n;
  }
  void push(Node n) { d_.push_back(n); }
  void clear() { d_.clear(); }

 private:
  std::vector<Node> d_;
};

template <>
class BasicQueue<true> {
 // Basically a queue
 public:
  BasicQueue() { begin_ = 0; }
  BasicQueue(IncBipart &) { begin_ = 0; }

  bool empty() const { return begin_ == d_.size(); }
  Node top() const { return d_[begin_]; }
  Node pop() {
    Node n = d_[begin_++];
    if (begin_ > d_.size() / 2) { // Half of the space wasted
      d_.erase(d_.begin(), d_.begin() + begin_);
      begin_ = 0;
    }
    return n;
  }
  void push(Node n) { d_.push_back(n); }
  void clear() { d_.clear(); begin_ = 0; }

 private:
  std::vector<Node> d_;
  std::size_t begin_;
};

// Queue that only handles positive and zero gain
template <bool FIFO=false>
class PosQueue {
 public:
  PosQueue(const IncBipart &inc)
    : inc_(inc) {}

  bool empty() {
    while (!pos_gain_.empty()) {
      Node n = pos_gain_.top();
      if (inc_.gain(n) <= 0) {
        pos_gain_.pop();
        zero_gain_.push(n);
      }
      else return false;
    }
    while (!zero_gain_.empty()) {
      Node n = zero_gain_.top();
      if (inc_.gain(n) < 0) {
        zero_gain_.pop();
      }
      else return false;
    }
    return true;
  }

  Node pop() {
    empty(); // Consume all out-of-place nodes
    if (!pos_gain_.empty()) {
      Node n = pos_gain_.top();
      pos_gain_.pop();
      return n;
    }
    Node n = zero_gain_.top();
    zero_gain_.pop();
    return n;
  }

  void push(Node n) {
    Weight g = inc_.gain(n);
    if (g > 0) {
      pos_gain_.push(n);
    }
    else if (g == 0) {
      zero_gain_.push(n);
    }
  }

  void clear() {
    pos_gain_.clear();
    zero_gain_.clear();
  }

 private:
  BasicQueue<FIFO> pos_gain_;
  BasicQueue<FIFO> zero_gain_;
  const IncBipart &inc_;
};

// Queue with positive, zero, and negative gain
template <bool FIFO=false>
class ThresholdQueue {
 public:
  ThresholdQueue(const IncBipart &inc)
    : inc_(inc) {}

  bool empty() {
    while (!pos_gain_.empty()) {
      Node n = pos_gain_.top();
      if (inc_.gain(n) <= 0) {
        pos_gain_.pop();
        zero_gain_.push(n);
      }
      else return false;
    }
    while (!zero_gain_.empty()) {
      Node n = zero_gain_.top();
      if (inc_.gain(n) < 0) {
        zero_gain_.pop();
        neg_gain_.push(n);
      }
      else return false;
    }
    return neg_gain_.empty();
  }

  Node pop() {
    while (!pos_gain_.empty()) {
      Node n = pos_gain_.top();
      pos_gain_.pop();
      if (inc_.gain(n) > 0) return n;
      else zero_gain_.push(n);
    }
    while (!zero_gain_.empty()) {
      Node n = zero_gain_.top();
      zero_gain_.pop();
      if (inc_.gain(n) >= 0) return n;
      else neg_gain_.push(n);
    }
    Node n = neg_gain_.top();
    neg_gain_.pop();
    return n;
  }

  void push(Node n) {
    Weight g = inc_.gain(n);
    if (g > 0) {
      pos_gain_.push(n);
    }
    else if (g == 0) {
      zero_gain_.push(n);
    }
    else {
      neg_gain_.push(n);
    }
  }

  void clear() {
    pos_gain_.clear();
    zero_gain_.clear();
    neg_gain_.clear();
  }

 private:
  BasicQueue<FIFO> pos_gain_;
  BasicQueue<FIFO> zero_gain_;
  BasicQueue<FIFO> neg_gain_;
  const IncBipart &inc_;
};

// Typical bucket queue
template <bool FIFO=false>
class FMQueue {
 public:
  FMQueue(const IncBipart &inc)
    : inc_(inc)
    , last_bucket_(0) {
    max_gain_ = 0;
    for (Node n : inc.nodes()) {
      Weight loc_gain = 0;
      for (Edge e : inc.edges(n)) {
        loc_gain += inc.hypergraph().weight(e);
      }
      max_gain_ = std::max(loc_gain, max_gain_);
    }
    buckets_.resize(2 * max_gain_ + 1);
  }

  bool empty() {
    while (true) {
      while (last_bucket_ != 0 && buckets_[last_bucket_].empty()) {
        --last_bucket_;
      }
      if (buckets_[last_bucket_].empty()) return true;
      Node n = buckets_[last_bucket_].top();
      Weight bucket = inc_.gain(n) + max_gain_;
      if (bucket != last_bucket_) {
        buckets_[last_bucket_].pop();
        buckets_[bucket].push(n);
        last_bucket_ = std::max(last_bucket_, bucket);
      }
      else return false;
    }
  }

  Node pop() {
    empty(); // Consume all out-of-place nodes
    Node n = buckets_[last_bucket_].top();
    buckets_[last_bucket_].pop();
    return n;
  }

  void push(Node n) {
    Weight bucket = max_gain_ + inc_.gain(n);
    buckets_[bucket].push(n);
    last_bucket_ = std::max(bucket, last_bucket_);
  }

  void clear() {
    for (auto & bucket : buckets_) {
      bucket.clear();
    }
  }

 private:
  const IncBipart &inc_;
  std::vector<BasicQueue<FIFO> > buckets_;
  Weight last_bucket_;
  Weight max_gain_;
};


} // End namespace minipart


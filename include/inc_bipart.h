// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#pragma once

#include "hypergraph.h"

#include <array>

namespace minipart {

template <class Index, class Weight, class Resource>
std::int64_t computeBipartCost(const Hypergraph<Index, Weight>&, const Mapping<std::uint8_t>&);

template <class Index, class Weight, class Resource>
class IncBipart {
  static_assert(std::is_signed<Weight>::value, "Gain type must be signed");
  static_assert(std::is_signed<Resource>::value, "Resource type must be signed");
  typedef std::array<Index, 2> CounterPair;

 public:
  IncBipart(const Problem<Index, Weight, Resource> &pb);

  std::size_t nNodes() const { return h_.nNodes(); }
  std::size_t nEdges() const { return h_.nEdges(); }
  std::size_t nResources() const { return capacities_.size2(); }

  Range<Node<Index> > nodes() const { return h_.nodes(); }

  Weight cost() const { return cost_; }
  bool legal() const;

  Weight gain(Node<Index> n) const { return gains_[n.id]; }
  bool mapping(Node<Index> n) const { return mapping_[n].id; }

  bool canMove(Node<Index> n) const;

  void move(Node<Index> n);
  template <typename F>
  void move(Node<Index> n, const F &onGainIncrease);

  bool tryMove(Node<Index> n);
  template <typename F>
  bool tryMove(Node<Index> n, const F &onGainIncrease);

  bool cut(Edge<Index> e) const { return edgeState_[e.id][0] != 0 && edgeState_[e.id][1] != 0; }
  bool overflow(bool partition) const;

  void checkConsistency() const;

 private:
  std::vector<CounterPair> initState() const;
  Weight initCost() const;
  std::vector<Weight> initGains() const;
  Matrix<Resource> initRemaining() const;

 private:
  Hypergraph<Index, Weight> h_;
  Mapping<std::uint8_t> mapping_;

  Matrix<Resource> demands_;
  Matrix<Resource> capacities_;
  Matrix<Resource> remaining_;

  std::vector<CounterPair> edgeState_;
  std::vector<Weight> gains_;
  Weight cost_;
};

template <class Index, class Weight, class Resource>
IncBipart<Index, Weight, Resource>::IncBipart(const Problem<Index, Weight, Resource> &pb)
: h_(pb.hypergraph)
, mapping_(pb.hypergraph.nNodes())
, demands_ (pb.demands)
, capacities_ (pb.capacities) {
  edgeState_ = initState();
  cost_ = initCost();
  gains_ = initGains();
  remaining_ = initRemaining();
}

template <class Index, class Weight, class Resource>
std::vector<typename IncBipart<Index, Weight, Resource>::CounterPair> IncBipart<Index, Weight, Resource>::initState() const {
  std::vector<CounterPair> ret(h_.nEdges(), CounterPair({0, 0}));
  for (auto e : h_.edges()) {
    Index count = 0;
    for (auto n : h_.nodes(e)) {
      count += mapping(n);
    }
    ret[e.id][0] = h_.nodes(e).size() - count;
    ret[e.id][1] = count;
  }
  return ret;
}

template <class Index, class Weight, class Resource>
bool IncBipart<Index, Weight, Resource>::legal() const {
  for (std::size_t i = 0; i < 2; ++i) {
    if (overflow(i)) return false;
  }
  return true;
}

template <class Index, class Weight, class Resource>
bool IncBipart<Index, Weight, Resource>::overflow(bool i) const {
  for (std::size_t j = 0; j < nResources(); ++j) {
    if (remaining_(i, j) < 0) return true;
  }
  return false;
}

template <class Index, class Weight, class Resource>
bool IncBipart<Index, Weight, Resource>::canMove(Node<Index> n) const {
  bool to = !mapping(n);
  for (std::size_t j = 0; j < nResources(); ++j) {
    if (remaining_(to, j) < demands_(n.id, j)) return false;
  }
  return true;
}

template <class Index, class Weight, class Resource>
Weight IncBipart<Index, Weight, Resource>::initCost() const {
  Weight ret = 0;
  for (auto e : h_.edges()) {
    CounterPair cnt = edgeState_[e.id];
    bool cut = (cnt[0] != 0) && (cnt[1] != 0);
    ret += h_.weight(e) * cut;
  }
  return ret;
}

template <class Index, class Weight, class Resource>
std::vector<Weight> IncBipart<Index, Weight, Resource>::initGains() const {
  std::vector<Weight> ret(h_.nNodes(), 0);
  for (auto n : h_.nodes()) {
    bool from = mapping(n);
    bool to = !from;
    Weight g = 0;
    for (auto e : h_.edges(n)) {
      CounterPair cnt = edgeState_[e.id];
      g += (cnt[from] == 1) * h_.weight(e);
      g -= (cnt[to]   == 0) * h_.weight(e);
    }
    ret[n.id] = g;
  }
  return ret;
}

template <class Index, class Weight, class Resource>
Matrix<Resource> IncBipart<Index, Weight, Resource>::initRemaining() const {
  Matrix<Resource> ret = capacities_;
  for (auto n : h_.nodes()) {
    bool part = mapping(n);
    for (std::size_t j = 0; j < nResources(); ++j) {
      ret(part, j) -= demands_(n.id, j);
    }
  }
  return ret;
}

template <class Index, class Weight, class Resource>
void IncBipart<Index, Weight, Resource>::move(Node<Index> n) {
  move(n, [](Node<Index> n, Weight w) {});
}

template <class Index, class Weight, class Resource>
template <typename F>
void IncBipart<Index, Weight, Resource>::move(Node<Index> n, const F &onGainIncrease) {
  bool from = mapping(n);
  bool to = !from;
  this->cost_ -= gains_[n.id];
  gains_[n.id] = -gains_[n.id];
  for (auto e : h_.edges(n)) {
    CounterPair &cnt = this->edgeState_[e.id];
    --cnt[from];
    ++cnt[to];

    if (cnt[from] == 0) {
      for (auto o : h_.nodes(e)) {
        if (o == n) continue;
        gains_[o.id] -= h_.weight(e);
      }
    }
    else if(cnt[from] == 1) {
      for (auto o : h_.nodes(e)) {
        if (o == n) continue;
        if (mapping(o) == to) continue;
        gains_[o.id] += h_.weight(e);
        onGainIncrease(o, gains_[o.id]);
      }
    }
    if (cnt[to] == 1) {
      for (auto o : h_.nodes(e)) {
        if (o == n) continue;
        gains_[o.id] += h_.weight(e);
        onGainIncrease(o, gains_[o.id]);
      }
    }
    else if (cnt[to] == 2) {
      for (auto o : h_.nodes(e)) {
        if (mapping(o) != to) continue;
        assert (o != n);
        gains_[o.id] -= h_.weight(e);
      }
    }
  }

  for (std::size_t j = 0; j < nResources(); ++j) {
    remaining_(to, j) -= demands_(n.id, j);
    remaining_(from, j) += demands_(n.id, j);
  }

  this->mapping_[n] = Part<std::uint8_t>(to);
}

template <class Index, class Weight, class Resource>
bool IncBipart<Index, Weight, Resource>::tryMove(Node<Index> n) {
  return tryMove(n, [](Node<Index> n, Weight w) {});
}

template <class Index, class Weight, class Resource>
template<typename F>
bool IncBipart<Index, Weight, Resource>::tryMove(Node<Index> n, const F &onGainIncrease) {
  if (canMove(n)) {
    move(n, onGainIncrease);
    return true;
  }
  return false;
}


template <class Index, class Weight, class Resource>
void IncBipart<Index, Weight, Resource>::checkConsistency() const {
  assert (mapping_.nNodes() == h_.nNodes());
  assert (demands_.size1() == nNodes());
  assert (capacities_.size1() == 2);
  assert (capacities_.size2() == demands_.size2());
  for (auto n : h_.nodes()) {
    assert (mapping_[n].id < 2);
  }
  assert (edgeState_ == initState());
  assert (cost_ == initCost());
  assert (initGains() == gains_);

  Matrix<Resource> check = initRemaining();
  for (std::size_t i = 0; i < 2; ++i) {
    for (std::size_t j = 0; j < nResources(); ++j) {
      assert (remaining_(i, j) == check(i, j));
    }
  }
}

template <class Index, class Weight, class Resource>
std::int64_t computeBipartCost(const Hypergraph<Index, Weight> &h, const Mapping<std::uint8_t>&m) {
  std::int64_t ret = 0;
  for (auto e : h.edges()) {
    // Used as a bitset
    std::uint8_t used = 0;
    for (auto n : h.nodes(e)) {
      uint8_t pos = m[n].id;
      assert (pos <= 1);
      used |= (pos + 1);
    }
    assert (used <= 3);
    // Cut if both bits are true
    if (used == 3) {
      ret += h.weight(e);
    }
  }
  return ret;
}


}  // End namespace minipart


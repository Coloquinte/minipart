// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#pragma once

#include "solver.h"

#include <array>
#include <random>

namespace minipart {

class IncBipart {
  static_assert(std::is_signed<Weight>::value, "Gain type must be signed");
  static_assert(std::is_signed<Resource>::value, "Resource type must be signed");
  typedef std::array<Index, 2> CounterPair;

 public:
  IncBipart(const Problem &pb);
  IncBipart(const Problem &pb, const Mapping &);

  std::size_t nNodes() const { return h_.nNodes(); }
  std::size_t nEdges() const { return h_.nEdges(); }
  std::size_t nResources() const { return capacities_.size2(); }

  Range<Node> nodes() const { return h_.nodes(); }
  Range<Edge> edges() const { return h_.edges(); }
  Slice<Node> nodes(Edge e) const { return h_.nodes(e); }
  Slice<Edge> edges(Node n) const { return h_.edges(n); }

  Weight cost() const { return cost_; }
  bool legal() const;
  Mapping mapping() const { return mapping_; }

  Weight gain(Node n) const { return gains_[n.id]; }
  bool mapping(Node n) const { return mapping_[n].id; }
  std::size_t moves() const { return moves_; } // Stats

  bool canMove(Node n) const;

  void move(Node n);
  template <typename F>
  void move(Node n, const F &onGainIncrease);

  bool tryMove(Node n);
  template <typename F>
  bool tryMove(Node n, const F &onGainIncrease);

  bool cut(Edge e) const;
  bool overflow(bool partition) const;

  const Hypergraph &hypergraph() const { return h_; }
  const Matrix<Resource> &demands() const { return demands_; }
  const Matrix<Resource> &capacities() const { return capacities_; }

  void checkConsistency() const;

 private:
  void init();
  std::vector<CounterPair> initState() const;
  Weight initCost() const;
  std::vector<Weight> initGains() const;
  Matrix<Resource> initRemaining() const;

 private:
  const Hypergraph &h_;
  Mapping mapping_;

  const Matrix<Resource> &demands_;
  const Matrix<Resource> &capacities_;
  Matrix<Resource> remaining_;

  std::vector<CounterPair> edgeState_;
  std::vector<Weight> gains_;
  Weight cost_;

  // Statistics
  std::size_t moves_;
};

// Basic high-level local search
void place(IncBipart &inc, std::minstd_rand &rgen, SolverOptions options);
void optimize(IncBipart &inc, std::minstd_rand &rgen, SolverOptions options);

inline IncBipart::IncBipart(const Problem &pb)
: h_(pb.hypergraph)
, mapping_(pb.hypergraph.nNodes())
, demands_ (pb.demands)
, capacities_ (pb.capacities) {
  init();
}

inline IncBipart::IncBipart(const Problem &pb, const Mapping &m)
: h_(pb.hypergraph)
, mapping_(m)
, demands_ (pb.demands)
, capacities_ (pb.capacities) {
  init();
}

inline void IncBipart::init() {
  edgeState_ = initState();
  cost_ = initCost();
  gains_ = initGains();
  remaining_ = initRemaining();
  moves_ = 0;
}

inline std::vector<typename IncBipart::CounterPair> IncBipart::initState() const {
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

inline bool IncBipart::legal() const {
  for (std::size_t i = 0; i < 2; ++i) {
    if (overflow(i)) return false;
  }
  return true;
}

inline bool IncBipart::cut(Edge e) const {
  return edgeState_[e.id][0] != 0
      && edgeState_[e.id][1] != 0;
}

inline bool IncBipart::overflow(bool i) const {
  for (std::size_t j = 0; j < nResources(); ++j) {
    if (remaining_(i, j) < 0) return true;
  }
  return false;
}

inline bool IncBipart::canMove(Node n) const {
  bool to = !mapping(n);
  for (std::size_t j = 0; j < nResources(); ++j) {
    if (remaining_(to, j) < demands_(n.id, j)) return false;
  }
  return true;
}

inline Weight IncBipart::initCost() const {
  Weight ret = 0;
  for (auto e : h_.edges()) {
    CounterPair cnt = edgeState_[e.id];
    bool cut = (cnt[0] != 0) && (cnt[1] != 0);
    ret += h_.weight(e) * cut;
  }
  return ret;
}

inline std::vector<Weight> IncBipart::initGains() const {
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

inline Matrix<Resource> IncBipart::initRemaining() const {
  Matrix<Resource> ret = capacities_;
  for (auto n : h_.nodes()) {
    bool part = mapping(n);
    for (std::size_t j = 0; j < nResources(); ++j) {
      ret(part, j) -= demands_(n.id, j);
    }
  }
  return ret;
}

inline void IncBipart::move(Node n) {
  move(n, [](Node n, Weight w) {});
}

template <typename F>
inline void IncBipart::move(Node n, const F &onGainIncrease) {
  ++moves_;
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

  this->mapping_[n] = Part(to);
}

inline bool IncBipart::tryMove(Node n) {
  return tryMove(n, [](Node n, Weight w) {});
}

template<typename F>
inline bool IncBipart::tryMove(Node n, const F &onGainIncrease) {
  if (canMove(n)) {
    move(n, onGainIncrease);
    return true;
  }
  return false;
}

inline void IncBipart::checkConsistency() const {
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

}  // End namespace minipart



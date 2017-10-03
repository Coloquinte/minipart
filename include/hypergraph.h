// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#pragma once

#include "common.h"

#include <unordered_set>

namespace minipart {

class HypergraphBuilder;

struct HEdgeData {
  Index limit_;
  Weight weight_;
  HEdgeData (Index l, Index w) : limit_(l), weight_(w) {}
};

class Hypergraph {
  static_assert(std::is_unsigned<Index>(), "Index type must be unsigned");
  static_assert(std::is_signed<Weight>(), "Weight type must be signed");

 public:
  typedef Node N;
  typedef Edge E;
  typedef HypergraphBuilder Builder;

  std::size_t nEdges() const;
  std::size_t nNodes() const;
  std::size_t nPins() const;

  Range<E> edges() const;
  Range<N> nodes() const;

  Slice<N> nodes   (Edge) const;
  Slice<E> edges   (Node) const;

  Weight weight(Edge) const;

  Hypergraph();
  Hypergraph(const Builder &);

  void checkConsistency() const;

 private:
  void checkLimits() const;
  void checkNodeToEdges() const;
  void checkUniquePins() const;

  void finalize();

 private:
  std::vector<HEdgeData > edges_;
  std::vector<N>    edge_pins_;

  std::vector<Index> node_limits_;
  std::vector<E> node_pins_;

  Index _nNodes;
};

class HypergraphBuilder {
 public:
  typedef Node N;
  typedef Edge E;
  typedef Hypergraph H;

  std::size_t nEdges() const;
  std::size_t nNodes() const;

  HypergraphBuilder(Index nNodes=0);
  N addNode();
  template <class It>
  E addEdge(It begin, It end, Weight w=1);
  E addEdge(std::initializer_list<N>, Weight w=1);

 private:
  std::vector<HEdgeData > edges_;
  std::vector<N>    edge_pins_;

  Index _nNodes;

  friend class Hypergraph;
};

struct Problem {
  Hypergraph hypergraph;
  Matrix<Resource> demands;
  Matrix<Resource> capacities;
};

inline Hypergraph::Hypergraph() {
  edges_.emplace_back(0, 0);
  _nNodes = 0;
}

inline Hypergraph::Hypergraph(const Builder &b) {
  edges_ = b.edges_;
  edge_pins_ = b.edge_pins_;
  _nNodes = b._nNodes;
  finalize();
}

inline std::size_t Hypergraph::nEdges() const {
  return edges_.size() - 1; // One dummy at the end for sparse storage consistency
}
inline std::size_t Hypergraph::nNodes() const {
  return _nNodes;
}
inline std::size_t Hypergraph::nPins() const {
  return edge_pins_.size();
}

inline Range<Edge > Hypergraph::edges() const {
  return Range<E>(E(0), E(nEdges()));
}
inline Range<Node > Hypergraph::nodes() const {
  return Range<N>(N(0), N(nNodes()));
}

inline Slice<Node > Hypergraph::nodes   (Edge e) const {
  return Slice<Node >(edge_pins_.begin() + edges_[e.id].limit_, edge_pins_.begin() + edges_[e.id+1].limit_);
}
inline Slice<Edge > Hypergraph::edges   (Node n) const {
  return Slice<Edge >(node_pins_.begin() + node_limits_[n.id], node_pins_.begin() + node_limits_[n.id+1]);
}

inline Weight Hypergraph::weight(Edge e) const {
  return edges_[e.id+1].weight_;
}

inline HypergraphBuilder::HypergraphBuilder(Index n) {
  edges_.emplace_back(0, 0);
  _nNodes = n;
}

inline std::size_t HypergraphBuilder::nEdges() const {
  return edges_.size() - 1; // One dummy at the end for sparse storage consistency
}

inline std::size_t HypergraphBuilder::nNodes() const {
  return _nNodes;
}

inline Node HypergraphBuilder::addNode() {
  return N(_nNodes++);
}

inline Edge HypergraphBuilder::addEdge(std::initializer_list<N> l, Weight w) {
  return addEdge(l.begin(), l.end(), w);
}

template<class It>
Edge HypergraphBuilder::addEdge(It begin, It end, Weight w) {
  E ret(nEdges());
  edge_pins_.insert(edge_pins_.end(), begin, end);
  edges_.emplace_back(edge_pins_.size(), w);
  return ret;
}

inline void Hypergraph::finalize() {
  edges_.shrink_to_fit();
  edge_pins_.shrink_to_fit();

  /*
   * Create the data for the nodes
   * We use _nodes for intermediate computations
   */

  node_limits_.assign(nNodes()+2, 0);
  node_pins_.resize(edge_pins_.size());

  // Get the number of Edges for each node
  for (E e : edges()) {
    for (N n : nodes(e)) {
      ++node_limits_[n.id+2];
    }
  }

  std::partial_sum(node_limits_.begin(), node_limits_.end(), node_limits_.begin());

  for (E e : edges()) {
    for (N n : nodes(e)) {
      node_pins_[node_limits_[n.id+1]++] = e;
    }
  }

  node_limits_.pop_back();
}

inline void Hypergraph::checkConsistency() const {
  checkLimits();
  checkNodeToEdges();
  checkUniquePins();
}

inline void Hypergraph::checkLimits() const {
  assert (nNodes() == node_limits_.size() - 1u);
  assert ((Index) nNodes() == _nNodes);
  assert (nEdges() ==  edges_.size() - 1u);
  assert (node_pins_.size() == edge_pins_.size());
  for (std::size_t i = 0; i < nNodes(); ++i) {
    assert (node_limits_[i] <= node_limits_[i+1]);
  }
  for (std::size_t i = 0; i < nEdges(); ++i) {
    assert (edges_[i].limit_ <= edges_[i+1].limit_);
  }

  assert (node_limits_.front() == (Index) 0);
  assert (node_limits_.back() == (Index) nPins());
}

inline void Hypergraph::checkNodeToEdges() const {
  std::vector<std::unordered_set<Index> > nodeToEdges(nNodes());

  for (N n : nodes()) {
    for (E e : edges(n)) {
      nodeToEdges[n.id].insert(e.id);
    }
  }

  for (E e : edges()) {
    for (N n : nodes(e)) {
      assert (nodeToEdges[n.id].count(e.id) != 0);
    }
  }
}

inline void Hypergraph::checkUniquePins() const {
  std::unordered_set<Index> pins;
  for (E e : edges()) {
    pins.clear();
    for (N n : nodes(e)) {
      assert (pins.count(n.id) == 0);
      pins.insert(n.id);
      assert (n.id >= 0);
      assert (n.id < (Index) nNodes());
    }
  }

  for (N n : nodes()) {
    pins.clear();
    for (E e : edges(n)) {
      assert (pins.count(e.id) == 0);
      pins.insert(e.id);
      assert (e.id >= 0);
      assert (e.id < (Index) nEdges());
    }
  }

}

}  // End namespace minipart


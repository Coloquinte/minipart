// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#pragma once

#include "common.h"

#include <unordered_set>

namespace minipart {

template <class Index, class Weight>
class HypergraphBuilder;

template <class Index, class Weight>
struct HEdgeData {
  Index limit_;
  Weight weight_;
  HEdgeData (Index l, Index w) : limit_(l), weight_(w) {}
};

template <class Index, class Weight>
class Hypergraph {
  static_assert(std::is_unsigned<Index>(), "Index type must be unsigned");
  static_assert(std::is_signed<Weight>(), "Weight type must be signed");

 public:
  typedef Node<Index> N;
  typedef Edge<Index> E;
  typedef HypergraphBuilder<Index, Weight> Builder;

  std::size_t nEdges() const;
  std::size_t nNodes() const;
  std::size_t nPins() const;

  Range<E> edges() const;
  Range<N> nodes() const;

  Slice<N> nodes   (Edge<Index>) const;
  Slice<E> edges   (Node<Index>) const;

  Weight weight(Edge<Index>) const;

  Hypergraph();
  Hypergraph(const Builder &);

  void checkConsistency() const;

 private:
  void checkLimits() const;
  void checkNodeToEdges() const;
  void checkUniquePins() const;

  void finalize();

 private:
  std::vector<HEdgeData<Index, Weight> > edges_;
  std::vector<N>    edge_pins_;

  std::vector<Index> node_limits_;
  std::vector<E> node_pins_;

  Index _nNodes;
};

template <class Index, class Weight>
class HypergraphBuilder {
 public:
  typedef Node<Index> N;
  typedef Edge<Index> E;
  typedef Hypergraph<Index, Weight> H;

  std::size_t nEdges() const;
  std::size_t nNodes() const;

  HypergraphBuilder(Index nNodes=0);
  N addNode();
  template <class It>
  E addEdge(It begin, It end, Weight w=1);
  E addEdge(std::initializer_list<N>, Weight w=1);

 private:
  std::vector<HEdgeData<Index, Weight> > edges_;
  std::vector<N>    edge_pins_;

  Index _nNodes;

  friend class Hypergraph<Index, Weight>;
};

template <class Index, class Weight>
inline Hypergraph<Index, Weight>::Hypergraph() {
  edges_.emplace_back(0, 0);
  _nNodes = 0;
}

template <class Index, class Weight>
inline Hypergraph<Index, Weight>::Hypergraph(const Builder &b) {
  edges_ = b.edges_;
  edge_pins_ = b.edge_pins_;
  _nNodes = b._nNodes;
  finalize();
}

template <class Index, class Weight>
inline std::size_t Hypergraph<Index, Weight>::nEdges() const {
  return edges_.size() - 1; // One dummy at the end for sparse storage consistency
}
template <class Index, class Weight>
inline std::size_t Hypergraph<Index, Weight>::nNodes() const {
  return _nNodes;
}
template <class Index, class Weight>
inline std::size_t Hypergraph<Index, Weight>::nPins() const {
  return edge_pins_.size();
}

template <class Index, class Weight>
inline Range<Edge<Index> > Hypergraph<Index, Weight>::edges() const {
  return Range<E>(E(0), E(nEdges()));
}
template <class Index, class Weight>
inline Range<Node<Index> > Hypergraph<Index, Weight>::nodes() const {
  return Range<N>(N(0), N(nNodes()));
}

template <class Index, class Weight>
inline Slice<Node<Index> > Hypergraph<Index, Weight>::nodes   (Edge<Index> e) const {
  return Slice<Node<Index> >(edge_pins_.begin() + edges_[e.id].limit_, edge_pins_.begin() + edges_[e.id+1].limit_);
}
template <class Index, class Weight>
inline Slice<Edge<Index> > Hypergraph<Index, Weight>::edges   (Node<Index> n) const {
  return Slice<Edge<Index> >(node_pins_.begin() + node_limits_[n.id], node_pins_.begin() + node_limits_[n.id+1]);
}

template <class Index, class Weight>
inline Weight Hypergraph<Index, Weight>::weight(Edge<Index> e) const {
  return edges_[e.id+1].weight_;
}

template <class Index, class Weight>
inline HypergraphBuilder<Index, Weight>::HypergraphBuilder(Index n) {
  edges_.emplace_back(0, 0);
  _nNodes = n;
}

template <class Index, class Weight>
inline std::size_t HypergraphBuilder<Index, Weight>::nEdges() const {
  return edges_.size() - 1; // One dummy at the end for sparse storage consistency
}

template <class Index, class Weight>
inline std::size_t HypergraphBuilder<Index, Weight>::nNodes() const {
  return _nNodes;
}

template <class Index, class Weight>
inline Node<Index> HypergraphBuilder<Index, Weight>::addNode() {
  return N(_nNodes++);
}

template <class Index, class Weight>
inline Edge<Index> HypergraphBuilder<Index, Weight>::addEdge(std::initializer_list<N> l, Weight w) {
  return addEdge(l.begin(), l.end(), w);
}

template <class Index, class Weight>
template<class It>
Edge<Index> HypergraphBuilder<Index, Weight>::addEdge(It begin, It end, Weight w) {
  E ret(nEdges());
  edge_pins_.insert(edge_pins_.end(), begin, end);
  edges_.emplace_back(edge_pins_.size(), w);
  return ret;
}

template <class Index, class Weight>
inline void Hypergraph<Index, Weight>::finalize() {
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

template <class Index, class Weight>
inline void Hypergraph<Index, Weight>::checkConsistency() const {
  checkLimits();
  checkNodeToEdges();
  checkUniquePins();
}

template <class Index, class Weight>
inline void Hypergraph<Index, Weight>::checkLimits() const {
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

template <class Index, class Weight>
inline void Hypergraph<Index, Weight>::checkNodeToEdges() const {
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

template <class Index, class Weight>
inline void Hypergraph<Index, Weight>::checkUniquePins() const {
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


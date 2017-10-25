// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#pragma once

#include "common.h"

namespace minipart {

class HypergraphBuilder;
class Hypergraph;

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
  std::size_t nPins() const;

  Range<E> edges() const;
  Slice<N> nodes   (Edge) const;
  Weight weight(Edge) const;

  HypergraphBuilder(Index nNodes=0);
  N addNode();
  template <class It>
  E addEdge(It begin, It end, Weight w=1);
  E addEdge(std::initializer_list<N>, Weight w=1);

  // Sort and remove duplicate pins
  void finalize();
  // Merge identical edges
  void vectorize();

 private:
  std::vector<HEdgeData > edges_;
  std::vector<N>    edge_pins_;

  Index _nNodes;

  friend class Hypergraph;
};

std::int64_t computeCostBipart(const Hypergraph&, const Mapping&);
std::vector<int> countCutsBipart(const Hypergraph &, const std::vector<Mapping> &);

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

inline std::size_t HypergraphBuilder::nPins() const {
  return edge_pins_.size();
}

inline Range<Edge > HypergraphBuilder::edges() const {
  return Range<E>(E(0), E(nEdges()));
}

inline Slice<Node > HypergraphBuilder::nodes(Edge e) const {
  return Slice<Node >(edge_pins_.begin() + edges_[e.id].limit_, edge_pins_.begin() + edges_[e.id+1].limit_);
}

inline Weight HypergraphBuilder::weight(Edge e) const {
  return edges_[e.id+1].weight_;
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

}  // End namespace minipart


// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#pragma once

#include "common.h"

namespace minipart {

class HypergraphBuilder;
class Hypergraph;

struct Edge2Data {
  Node source;
  Node target;
  Weight weight;
  Edge2Data(Node s, Node t, Weight w) : source(s), target(t), weight(w) {}
};

struct Pin2Data {
  Node target;
  Weight weight;
};

class Hypergraph {
  static_assert(std::is_unsigned<Index>(), "Index type must be unsigned");
  static_assert(std::is_signed<Weight>(), "Weight type must be signed");

  /*
   * Edges with only 2 pins are treated entirely differently
   */

  struct EdgeData {
    Index limit_;
    Weight weight_;
    EdgeData (Index l, Weight w) : limit_(l), weight_(w) {}
  };

 public:
  typedef HypergraphBuilder Builder;

  std::size_t nNodes() const;
  std::size_t nEdges() const;

  Range<Node>  nodes() const;
  Range<Edge>  edges() const;
  const std::vector<Edge2Data> &edges2() const;

  Slice<Node>  nodes  (Edge) const;
  Slice<Edge>  edges  (Node) const;
  Weight weight(Edge) const;

  Slice<Pin2Data> edges2 (Node) const;

  Hypergraph();
  Hypergraph(const Builder &);

  void checkConsistency() const;

 private:
  void checkLimits() const;
  void checkNodeToEdges() const;
  void checkUniquePins() const;

  void finalize();

 private:
  std::vector<EdgeData > edges_;
  std::vector<Node>    edge_pins_;

  std::vector<Index> node_limits_;
  std::vector<Edge> node_pins_;

  std::vector<Index> node_limits2_;
  std::vector<Pin2Data> node_pins2_;
  std::vector<Edge2Data> edge_pins2_;

  Index _nNodes;

  friend class HypergraphBuilder;
};

class HypergraphBuilder {
 public:
  typedef Hypergraph H;

  std::size_t nEdges() const;
  std::size_t nNodes() const;

  Range<Edge> edges() const;
  Slice<Node> nodes   (Edge) const;
  Weight weight(Edge) const;

  HypergraphBuilder(Index nNodes=0);
  Node addNode();
  template <class It>
  Edge addEdge(It begin, It end, Weight w=1);
  Edge addEdge(std::initializer_list<Node>, Weight w=1);
  Edge2 addEdge2(Node a, Node b, Weight w=1);

  // Sort and remove duplicate pins
  void finalize();
  // Merge identical edges
  void vectorize();

 private:
  std::vector<H::EdgeData> edges_;
  std::vector<Node>    edge_pins_;
  std::vector<Edge2Data> edge_pins2_;

  Index _nNodes;

  friend class Hypergraph;
};

std::int64_t computeCostBipart(const Hypergraph&, const Mapping&);
std::int64_t computeCostCut(const Hypergraph&, const Mapping&);
std::int64_t computeCostDegree(const Hypergraph&, const Mapping&);

std::vector<int> countCutsBipart(const Hypergraph &, const std::vector<Mapping> &);
std::vector<int> countCuts(const Hypergraph &, const std::vector<Mapping> &);

inline Hypergraph::Hypergraph() {
  edges_.emplace_back(0, 0);
  _nNodes = 0;
}

inline Hypergraph::Hypergraph(const Builder &b) {
  edges_ = b.edges_;
  edge_pins_ = b.edge_pins_;
  edge_pins2_ = b.edge_pins2_;
  _nNodes = b._nNodes;
  finalize();
}

inline std::size_t Hypergraph::nEdges() const {
  return edges_.size() - 1; // One dummy at the end for sparse storage consistency
}
inline std::size_t Hypergraph::nNodes() const {
  return _nNodes;
}

inline Range<Edge > Hypergraph::edges() const {
  return Range<Edge>(Edge(0), Edge(nEdges()));
}
inline Range<Node > Hypergraph::nodes() const {
  return Range<Node>(Node(0), Node(nNodes()));
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

inline Range<Edge > HypergraphBuilder::edges() const {
  return Range<Edge>(Edge(0), Edge(nEdges()));
}

inline Slice<Node > HypergraphBuilder::nodes(Edge e) const {
  return Slice<Node >(edge_pins_.begin() + edges_[e.id].limit_, edge_pins_.begin() + edges_[e.id+1].limit_);
}

inline Weight HypergraphBuilder::weight(Edge e) const {
  return edges_[e.id+1].weight_;
}

inline Node HypergraphBuilder::addNode() {
  return Node(_nNodes++);
}

inline Edge HypergraphBuilder::addEdge(std::initializer_list<Node> l, Weight w) {
  return addEdge(l.begin(), l.end(), w);
}

template<class It>
inline Edge HypergraphBuilder::addEdge(It begin, It end, Weight w) {
  Edge ret(nEdges());
  edge_pins_.insert(edge_pins_.end(), begin, end);
  edges_.emplace_back(edge_pins_.size(), w);
  return ret;
}

inline Edge2 HypergraphBuilder::addEdge2(Node s, Node t, Weight w) {
  Edge2 ret(edge_pins2_.size());
  Index mn = std::min(s.id, t.id);
  Index mx = std::max(s.id, t.id);
  edge_pins2_.emplace_back(Node(mn), Node(mx), w);
  return ret;
}

inline Slice<Pin2Data> Hypergraph::edges2(Node n) const {
  return Slice<Pin2Data>(node_pins2_.begin() + node_limits2_[n.id], node_pins2_.begin() + node_limits2_[n.id+1]);
}

inline const std::vector<Edge2Data> &Hypergraph::edges2() const {
  return edge_pins2_;
}

}  // End namespace minipart


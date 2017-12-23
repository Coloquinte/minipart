// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#include "hypergraph.h"

#include <unordered_set>
#include <algorithm>

namespace minipart {

void Hypergraph::finalize() {
  // Finalize large edges
  edges_.shrink_to_fit();
  edge_pins_.shrink_to_fit();

  /*
   * Create the data for the nodes
   * We use nodes_limits_ for intermediate computations
   */

  node_limits_.assign(nNodes()+2, 0);
  node_pins_.resize(edge_pins_.size());

  // Get the number of Edges for each node
  for (Edge e : edges()) {
    for (Node n : nodes(e)) {
      assert (n.id < nNodes());
      ++node_limits_[n.id+2];
    }
  }

  std::partial_sum(node_limits_.begin(), node_limits_.end(), node_limits_.begin());

  for (Edge e : edges()) {
    for (Node n : nodes(e)) {
      node_pins_[node_limits_[n.id+1]++] = e;
    }
  }
  node_limits_.pop_back();

  // Finalize 2-Edges
  node_limits2_.assign(nNodes()+2, 0);
  node_pins2_.resize(2 * edge_pins2_.size());
  for (auto e : edge_pins2_) {
    ++node_limits2_[e.source.id+2];
    ++node_limits2_[e.target.id+2];
  }

  std::partial_sum(node_limits2_.begin(), node_limits2_.end(), node_limits2_.begin());

  for (auto e : edge_pins2_) {
      node_pins2_[node_limits2_[e.source.id+1]++] = {e.target, e.weight};
      node_pins2_[node_limits2_[e.target.id+1]++] = {e.source, e.weight};
  }
  node_limits2_.pop_back();
}

void Hypergraph::checkConsistency() const {
  checkLimits();
  checkNodeToEdges();
  checkUniquePins();
}

void Hypergraph::checkLimits() const {
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
  assert (node_limits_.back() == (Index) edge_pins_.size());
}

void Hypergraph::checkNodeToEdges() const {
  std::vector<std::unordered_set<Index> > nodeToEdges(nNodes());

  for (Node n : nodes()) {
    for (Edge e : edges(n)) {
      nodeToEdges[n.id].insert(e.id);
    }
  }

  for (Edge e : edges()) {
    for (Node n : nodes(e)) {
      assert (nodeToEdges[n.id].count(e.id) != 0);
    }
  }
}

void Hypergraph::checkUniquePins() const {
  std::unordered_set<Index> pins;
  for (Edge e : edges()) {
    pins.clear();
    for (Node n : nodes(e)) {
      assert (pins.count(n.id) == 0);
      pins.insert(n.id);
      assert (n.id < (Index) nNodes());
    }
  }

  for (auto e : edges2()) {
    assert (e.source.id < (Index) nNodes());
    assert (e.target.id < (Index) nNodes());
  }

  for (Node n : nodes()) {
    pins.clear();
    for (Edge e : edges(n)) {
      assert (pins.count(e.id) == 0);
      pins.insert(e.id);
      assert (e.id >= 0);
      assert (e.id < (Index) nEdges());
    }
  }
}

void HypergraphBuilder::finalize() {
  // TODO: in-place?
  HypergraphBuilder ret(nNodes());
  for (Edge edge : edges()) {
    auto b = ret.edge_pins_.insert(ret.edge_pins_.end(), nodes(edge).begin(), nodes(edge).end());
    auto e = ret.edge_pins_.end();
    // Sort-uniq
    std::sort(b, e);
    auto it = std::unique(b, e);
    // Remove empty edge
    if (it - b <= 1) {
      ret.edge_pins_.resize(b - ret.edge_pins_.begin());
    }
    else {
      ret.edge_pins_.resize(it - ret.edge_pins_.begin());
      ret.edges_.emplace_back(ret.edge_pins_.size(), weight(edge));
    }
  }
  ret.edge_pins2_ = edge_pins2_;
  std::swap(*this, ret);
}

void HypergraphBuilder::vectorize() {
  // TODO: in-place?
  finalize();

  // Likely faster to separate the processing by size
  std::vector<std::vector<Edge> > size_to_edges(40);
  for (Edge e : edges()) {
    Index edge_size = nodes(e).size();
    if (edge_size >= size_to_edges.size()) size_to_edges.resize(edge_size+1);
    size_to_edges[edge_size].push_back(e);
    assert (std::is_sorted(nodes(e).begin(), nodes(e).end()));
  }

  // Merge identical edges
  HypergraphBuilder ret(nNodes());
  for (Index sz = 3; sz < size_to_edges.size(); ++sz) {
    const std::vector<Edge> &v = size_to_edges[sz];
    if (v.empty()) continue;

    auto edges_less = [&](Edge a, Edge b) -> bool {
      auto b1 = nodes(a).begin();
      auto b2 = nodes(b).begin();
      return std::lexicographical_compare(b1, b1+sz, b2, b2+sz);
    };
    std::sort(size_to_edges[sz].begin(), size_to_edges[sz].end(), edges_less);

    auto edges_equal = [&](Edge a, Edge b) -> bool {
      auto b1 = nodes(a).begin();
      auto b2 = nodes(b).begin();
      return std::equal(b1, b1+sz, b2, b2+sz);
    };

    Edge cur = v[0];
    Weight cur_weight = weight(cur);
    for (Index j = 1; j < v.size(); ++j) {
      if (!edges_equal(cur, v[j])) {
        auto b = nodes(cur).begin();
        ret.addEdge(b, b+sz, cur_weight);
        cur = v[j];
        cur_weight = weight(cur);
      }
      else {
        cur_weight += weight(v[j]);
      }
    }
    auto b = nodes(cur).begin();
    ret.addEdge(b, b+sz, cur_weight);
  }

  for (Edge e : size_to_edges[2]) {
    addEdge2(nodes(e)[0], nodes(e)[1], weight(e));
  }
  auto edges2_less = [](Edge2Data a, Edge2Data b) -> bool {
    return a.source < b.source
        || (a.source == b.source && a.target < b.target);
  };
  std::sort(edge_pins2_.begin(), edge_pins2_.end(), edges2_less);
  auto edges2_equal = [](Edge2Data a, Edge2Data b) -> bool {
    return a.source == b.source && a.target == b.target;
  };

  if (!edge_pins2_.empty()) {
    Edge2Data cur = edge_pins2_[0];
    for (Index j = 1; j < edge_pins2_.size(); ++j) {
      if (!edges2_equal(cur, edge_pins2_[j])) {
        ret.edge_pins2_.push_back(cur);
        cur = edge_pins2_[j];
      }
      else {
        cur.weight += edge_pins2_[j].weight;
      }
    }
    ret.edge_pins2_.push_back(cur);
  }

  std::swap(*this, ret);
}

}  // End namespace minipart



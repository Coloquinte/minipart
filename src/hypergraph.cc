// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#include "hypergraph.h"

#include <unordered_set>
#include <algorithm>

namespace minipart {

void Hypergraph::finalize() {
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
      assert (n.id < nNodes());
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
  assert (node_limits_.back() == (Index) nPins());
}

void Hypergraph::checkNodeToEdges() const {
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

void Hypergraph::checkUniquePins() const {
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

  // Get identical edges together
  for (Index sz = 2; sz < size_to_edges.size(); ++sz) {
    auto compare_edges = [&](Edge a, Edge b) -> bool {
      auto b1 = nodes(a).begin();
      auto b2 = nodes(b).begin();
      return std::lexicographical_compare(b1, b1+sz, b2, b2+sz);
    };
    std::sort(size_to_edges[sz].begin(), size_to_edges[sz].end(), compare_edges);
  }

  // Merge identical edges
  HypergraphBuilder ret(nNodes());
  for (Index sz = 2; sz < size_to_edges.size(); ++sz) {
    const std::vector<Edge> &v = size_to_edges[sz];
    if (v.empty()) continue;

    auto compare_edges = [&](Edge a, Edge b) -> bool {
      auto b1 = nodes(a).begin();
      auto b2 = nodes(b).begin();
      return std::equal(b1, b1+sz, b2, b2+sz);
    };
    Edge cur = v[0];
    Weight cur_weight = weight(cur);
    for (Index j = 1; j < v.size(); ++j) {
      if (!compare_edges(cur, v[j])) {
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

  std::swap(*this, ret);
}

}  // End namespace minipart



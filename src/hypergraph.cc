// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#include "hypergraph.h"

#include <unordered_set>

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

}  // End namespace minipart



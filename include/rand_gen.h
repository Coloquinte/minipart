// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#pragma once

#include "hypergraph.h"

#include <random>

namespace minipart {

template<typename Index, typename Weight, typename Rgen>
void addRandomEdges(HypergraphBuilder<Index, Weight> &h, Index nEdges, double average_degree, double average_weight, Rgen &rgen, Index min_degree = 2) {
  if (h.nNodes() == 0) return;

  assert (average_degree > min_degree);

  std::geometric_distribution<Index> degree_dist(1.0 / (1.0 + average_degree - min_degree) );
  std::geometric_distribution<Weight> weight_dist(1.0 / average_weight);
  std::uniform_int_distribution<Index> node_dist(0, h.nNodes()-1);

  std::unordered_set<Index> pinSet;
  std::vector<Node<Index> > pins;
  for (Index i = 0; i < nEdges; ++i) {
    pinSet.clear();
    pins.clear();
    Index sz = degree_dist(rgen) + min_degree;
    Weight w = 1 + weight_dist(rgen);
    for (Index j = 0; j < sz; ++j) {
      pinSet.insert(node_dist(rgen));
    }
    if (pinSet.size() < min_degree) continue;
    for (Index p : pinSet) {
      pins.emplace_back(p);
    }
    h.addEdge(pins.begin(), pins.end(), w);
  }
}

}  // End namespace minipart



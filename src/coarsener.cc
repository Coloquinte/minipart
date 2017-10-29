// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#include "coarsener.h"

#include <unordered_map>

namespace minipart {

Mapping Coarsening::operator() (const Mapping &m) const {
  Mapping ret(nNodesOut());
  assert (nNodesIn() == m.nNodes());
  for (std::size_t i = 0; i < nNodesIn(); ++i) {
    Node n(i);
    Node c = m_[i];
    ret[c] = m[n];
  }

  // Check that the coarsening is consistent with the mapping
  // It should only merge nodes that belong to the same partition
  for (std::size_t i = 0; i < nNodesIn(); ++i) {
    Node n(i);
    Node c = m_[i];
    assert(ret[c] == m[n]);
  }
  return ret;
}

Problem Coarsening::operator() (const Problem &pb) const {
  Problem ret;

  // Coarsen the demands
  ret.demands = Matrix<Resource>(nNodesOut(), pb.demands.size2());
  for (std::size_t i = 0; i < nNodesOut(); ++i) {
    for (std::size_t j = 0; j < pb.demands.size2(); ++j) {
      ret.demands(i, j) = 0;
    }
  }
  for (std::size_t i = 0; i < nNodesIn(); ++i) {
    for (std::size_t j = 0; j < pb.demands.size2(); ++j) {
      Node c = m_[i];
      ret.demands(c.id, j) += pb.demands(i, j);
    }
  }

  // Capacities are the same
  ret.capacities = pb.capacities;

  // Coarsen and simplify the hypergraph
  HypergraphBuilder b(nNodesOut());
  std::vector<Node> nodes;
  for (Edge e : pb.hypergraph.edges()) {
    for (Node n : pb.hypergraph.nodes(e)) {
      nodes.push_back(m_[n.id]);
    }
    b.addEdge(nodes.begin(), nodes.end(), pb.hypergraph.weight(e));
    nodes.clear();
  }
  b.vectorize();
  ret.hypergraph = b;

  return ret;
}

Mapping Coarsening::reverse (const Mapping &m) const {
  Mapping ret(nNodesIn());
  assert (nNodesOut() == m.nNodes());
  for (std::size_t i = 0; i < nNodesIn(); ++i) {
    Node n(i);
    Node c = m_[i];
    ret[n] = m[c];
  }
  return ret;
}

Coarsening infer_coarsening_blackbox(const std::vector<Mapping> &mappings) {
  // Coarsen nodes that are always together in the same partition
  assert (mappings.size() <= 64);

  std::vector<std::uint64_t> placements(mappings.front().nNodes(), 0);
  for (std::size_t i = 0; i < placements.size(); ++i) {
    std::uint64_t bs = 0;
    for (const Mapping &m : mappings) {
      bs = (bs << 1) | m[Node(i)].id;
    }
    placements[i] = bs;
  }

  std::size_t cur_id = 0;
  std::unordered_map<std::uint64_t, std::size_t> placement_to_coarsening;
  placement_to_coarsening.reserve(placements.size());
  for (auto pl : placements) {
    bool inserted = placement_to_coarsening.emplace(pl, cur_id).second;
    if (inserted) ++cur_id;
  }
  Coarsening ret(placements.size(), cur_id+1);
  for (std::size_t i = 0; i < placements.size(); ++i) {
    std::size_t id = placement_to_coarsening.at(placements[i]);
    ret[Node(i)] = Node(id);
  }
  return ret;
}

Coarsening infer_coarsening_connected(const Hypergraph &h, const std::vector<Mapping> &mappings) {
  // Coarsen edges that are never cut
  std::vector<int> cuts = countCutsBipart(h, mappings);

  //std::vector<Index> ids(h.nNodes(), 0);
  //Index cur_id = 0;
  Coarsening ret(h.nNodes(), 0);

  // Typical DFS for connected components
  std::vector<Node> nodes_to_visit;
  std::vector<char> node_visited(h.nNodes(), 0);
  std::vector<char> edge_visited(h.nEdges(), 0);

  for (Node n : h.nodes()) {
    if (node_visited[n.id]) continue;

    Node o = ret.addNode();
    nodes_to_visit.push_back(n);

    while (!nodes_to_visit.empty()) {
      Node v = nodes_to_visit.back();
      nodes_to_visit.pop_back();

      if (node_visited[v.id]) continue;
      node_visited[v.id] = true;
      ret[v] = o;

      for (Edge e : h.edges(v)) {
        // This avoids quadratic complexity
        if (edge_visited[e.id]) continue;
        edge_visited[e.id] = true;

        // Only coarsen edges which are never cut
        if (cuts[e.id] != 0) continue;

        for (Node s : h.nodes(e)) {
          nodes_to_visit.push_back(s);
        }
      }
    }

    assert (nodes_to_visit.empty());
  }

  return ret;
}

std::pair<Coarsening, std::vector<Mapping> > select_for_coarsening(const std::vector<Mapping> &mappings, std::size_t target_nnodes) {
  // Select only the first few mappings to get a coarse solution
  //   i.e. the largest solution that has that many nodes or fewer
  // Assume the mapping are already sorted in a nice order - best first for example
  assert (mappings.size() <= 64);

  std::vector<Mapping> selected(mappings.begin(), mappings.end());
  Coarsening ret;

  for (std::size_t i = 1; i < mappings.size(); ++i) {
    ret = infer_coarsening_blackbox(selected);
    if (ret.nNodesOut() <= target_nnodes) break;
    selected.pop_back();
  }
  return std::make_pair(ret, selected);
}

// Give several coarsening solutions, each associated with part of the pool
std::vector<std::pair<Coarsening, std::vector<Mapping> > > select_pool_coarsenings(const Problem &pb, const std::vector<Mapping> &pool, std::size_t target_n_nodes, std::minstd_rand &rgen) {
  std::vector<std::pair<Coarsening, std::vector<Mapping> > > ret;

  // Try using all solutions
  Coarsening using_all = infer_coarsening_blackbox(pool);
  if (using_all.nNodesOut() < 0.8 * pb.hypergraph.nNodes()) {
    // If the coarsening is good
    ret.emplace_back(using_all, pool);
    return ret;
  }

  // The problem size was not reduced enough
  // Use two coarsenings with half the pool each
  std::vector<Mapping> shuffled = pool;
  std::shuffle(shuffled.begin(), shuffled.end(), rgen);
  for (std::size_t i = 0; i < 2; ++i) {
    std::vector<Mapping> selected;
    for (std::size_t j = i; j < shuffled.size(); j += 2) {
      selected.push_back(shuffled[j]);
    }
    ret.emplace_back(infer_coarsening_blackbox(selected), selected);
  }

  return ret;
}

}  // End namespace minipart


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
  HypergraphBuilder hb(nNodesOut());
  std::vector<Node> nodes;
  for (Edge e : pb.hypergraph.edges()) {
    for (Node n : pb.hypergraph.nodes(e)) {
      nodes.push_back(m_[n.id]);
    }
    hb.addEdge(nodes.begin(), nodes.end(), pb.hypergraph.weight(e));
    nodes.clear();
  }
  for (auto e : pb.hypergraph.edges2()) {
    // TODO: use addEdge2 directly and include it in the vectorization
    hb.addEdge({m_[e.source.id], m_[e.target.id]}, e.weight);
  }
  hb.vectorize();

  ret.hypergraph = hb;
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

void Coarsening::checkConsistency() const {
  std::vector<char> used(nNodesOut(), 0);
  for (Node n : m_) {
    assert (n.id < nNodesOut());
    used[n.id] = true;
  }
  for (char u : used) {
    assert (u);
  }
}

bool mappings_are_bipart(const std::vector<Mapping> &mappings) {
  bool ret = true;
  for (const Mapping &m : mappings) {
    for (Index i = 0; i < m.nNodes(); ++i) {
      ret |= (m[Node(i)].id > 1);
    }
  }
  return ret;
}

Coarsening infer_coarsening_blackbox(const Problem &pb, const std::vector<Mapping> &mappings) {
  // Coarsen nodes that are always together in the same partition
  // TODO: make it faster and generalize to non-bipartitioning
  assert (mappings_are_bipart(mappings));
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
  Coarsening ret(placements.size(), cur_id);
  for (std::size_t i = 0; i < placements.size(); ++i) {
    std::size_t id = placement_to_coarsening.at(placements[i]);
    ret[Node(i)] = Node(id);
  }
  ret.checkConsistency();
  return ret;
}

Coarsening infer_coarsening_connected(const Problem &pb, const std::vector<Mapping> &mappings) {
  const Hypergraph &h = pb.hypergraph;

  // Coarsen edges that are never cut
  std::vector<int> cuts = countCutsBipart(h, mappings);
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

      // Push all neighbours on the stack
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

      // Same for 2-edges
      for (auto e : h.edges2(v)) {
        Node s = e.target;
        bool cut = false;
        for (const Mapping &m : mappings) {
          if (m[s] != m[v]) cut = true;
        }
        if (cut) continue;
        nodes_to_visit.push_back(s);
      }
    }

    assert (nodes_to_visit.empty());
  }

  ret.checkConsistency();
  return ret;
}

std::pair<Coarsening, std::vector<Mapping> > select_for_coarsening(const Problem &pb, const std::vector<Mapping> &mappings, std::size_t target_nnodes) {
  // Select only the first few mappings to get a coarse solution
  //   i.e. the largest solution that has that many nodes or fewer
  // Assume the mapping are already sorted in a nice order - best first for example
  assert (mappings.size() <= 64);

  std::vector<Mapping> selected(mappings.begin(), mappings.end());
  Coarsening ret;

  for (std::size_t i = 1; i < mappings.size(); ++i) {
    ret = infer_coarsening_blackbox(pb, selected);
    if (ret.nNodesOut() <= target_nnodes) break;
    selected.pop_back();
  }
  return std::make_pair(ret, selected);
}

// Give several coarsening solutions, each associated with part of the pool
// Criterions are:
//    * Less nodes --> smaller problems
//    * Sharing --> less work to do; more information on what's worth coarsening

std::vector<std::pair<Coarsening, std::vector<Mapping> > > select_pool_coarsenings(const Problem &pb, const std::vector<Mapping> &pool, std::size_t target_n_nodes, std::minstd_rand &rgen) {
  std::vector<std::pair<Coarsening, std::vector<Mapping> > > ret;

  // Try several techniques until we have something good enough
  Coarsening using_all = infer_coarsening_connected(pb, pool);
  if (using_all.nNodesOut() <= target_n_nodes) {
    ret.emplace_back(using_all, pool);
    return ret;
  }

  // Pure blackbox
  if (pool.size() <= 64) {
    Coarsening blackbox = infer_coarsening_blackbox(pb, pool);
    //assert (blackbox.nNodesOut() <= using_all.nNodesOut());
    if (blackbox.nNodesOut() <= target_n_nodes || pool.size() == 1) {
      ret.emplace_back(blackbox, pool);
      return ret;
    }
  }

  // If neither is coarse enough use two coarsenings with half the pool each
  std::vector<Mapping> shuffled = pool;
  std::shuffle(shuffled.begin(), shuffled.end(), rgen);
  for (std::size_t i = 0; i < 2; ++i) {
    std::vector<Mapping> selected;
    for (std::size_t j = i; j < shuffled.size(); j += 2) {
      selected.push_back(shuffled[j]);
    }
    ret.emplace_back(infer_coarsening_connected(pb, selected), selected);
  }

  return ret;
}

}  // End namespace minipart


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
  // TODO: non-blackbox algorithm, restricting coarse nodes to connected components
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

std::vector<std::pair<Coarsening, std::vector<Mapping> > > select_pool_coarsenings(const Problem &pb, const std::vector<Mapping> &pool, std::size_t target_n_nodes, std::minstd_rand &rgen) {
  std::vector<std::pair<Coarsening, std::vector<Mapping> > > ret;

  std::vector<Mapping> selected = pool;
  std::shuffle(selected.begin() + 1, selected.end(), rgen);
  ret.push_back(select_for_coarsening(selected, target_n_nodes));

  return ret;
}

}  // End namespace minipart


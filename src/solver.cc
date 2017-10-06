// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#include "inc_bipart.h"
#include "coarsener.h"

namespace minipart {

void place(IncBipart &inc, std::minstd_rand &rgen) {
  std::bernoulli_distribution dist;
  for (auto n : inc.nodes()) {
    if (dist(rgen)) inc.move(n);
  }

  std::vector<Node> nodes (inc.nodes().begin(), inc.nodes().end());
  std::shuffle(nodes.begin(), nodes.end(), rgen);
  for (auto n : nodes) {
    bool m = inc.mapping(n);
    if (inc.overflow(m)) inc.move(n);
  }

  inc.checkConsistency();
}

void optimize(IncBipart &inc, std::minstd_rand &rgen) {
  std::vector<Node > nodes (inc.nodes().begin(), inc.nodes().end());
  std::shuffle(nodes.begin(), nodes.end(), rgen);

  std::vector<Node > active_set, zero_gain_set;
  std::vector<int> num_zero_gain_moves(inc.nNodes(), 0);
  const int max_zero_gain_moves = 2;

  auto tryPush = [&](Node o) {
    if (inc.gain(o) > 0) {
      active_set.push_back(o);
    }
    else if (inc.gain(o) == 0 && ++num_zero_gain_moves[o.id] <= max_zero_gain_moves) {
      zero_gain_set.push_back(o);
    }
  };

  for (auto n : nodes) {
    tryPush(n);
  }

  while (!active_set.empty() || !zero_gain_set.empty()) {
    std::vector<Node > &pick = active_set.empty() ? zero_gain_set : active_set;
    Node n = pick.back();
    pick.pop_back();
    if (inc.gain(n) < 0) continue;
    inc.tryMove(n, [&](Node o, int w) { tryPush(o); });
  }
}

void coarsen_recurse(const Problem &pb, std::vector<Mapping> &mappings, std::minstd_rand &rgen);
void optimize_recurse(const Problem &pb, std::vector<Mapping> &mappings, std::minstd_rand &rgen);

void coarsen_recurse(const Problem &pb, std::vector<Mapping> &mappings, std::minstd_rand &rgen) {
  Coarsening coarsening = inferCoarsening(mappings);
  if (coarsening.nNodesOut() < 10 || coarsening.nNodesOut() >= 0.95 * pb.hypergraph.nNodes()) {
    return;
  }

  Problem c_pb = coarsening(pb);
  std::vector<Mapping> c_mappings;
  
  for (const Mapping &m : mappings) {
    auto c_m = coarsening(m);
    assert (computeBipartCost(pb.hypergraph, m) == computeBipartCost(c_pb.hypergraph, c_m));
    c_mappings.push_back(c_m);
  }
  optimize_recurse(c_pb, c_mappings, rgen);
  for (std::size_t i = 0; i < mappings.size(); ++i) {
    mappings[i] = coarsening.reverse(c_mappings[i]);
  }
}

void optimize_recurse(const Problem &pb, std::vector<Mapping> &mappings, std::minstd_rand &rgen) {
  for (Mapping &m : mappings) {
    IncBipart inc(pb, m);
    optimize(inc, rgen);
    m = inc.mapping();
  }
  coarsen_recurse(pb, mappings, rgen);
}

std::vector<Mapping> solve(const Problem &pb, int n_starts) {
  std::minstd_rand rgen;
  IncBipart inc(pb);

  std::vector<Mapping> mappings;

  for (int i = 0; i < n_starts; ++i) {
    place(inc, rgen);
    if (!inc.legal()) continue;
    optimize(inc, rgen);
    mappings.push_back(inc.mapping());
  }

  coarsen_recurse(pb, mappings, rgen);
  const int n_cycles = 3;
  for (int i = 0; i < n_cycles; ++i) {
    optimize_recurse(pb, mappings, rgen);
  }
  return mappings;
}

} // End namespace minipart


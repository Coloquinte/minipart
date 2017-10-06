// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#include "inc_bipart.h"

#include <iostream>
#include <iomanip>

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

  return mappings;
}

} // End namespace minipart


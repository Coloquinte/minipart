// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#include "inc_bipart.h"
#include "coarsener.h"

namespace minipart {

class ThresholdQueue {
 public:
  ThresholdQueue(const IncBipart &inc)
    : inc_(inc) {}

  bool empty() const {
    return pos_gain_.empty()
        && zero_gain_.empty()
        && neg_gain_.empty();
  }

  Node pop() {
    while (!pos_gain_.empty()) {
      Node n = pos_gain_.back();
      pos_gain_.pop_back();
      if (inc_.gain(n) > 0) return n;
      else zero_gain_.push_back(n);
    }
    while (!zero_gain_.empty()) {
      Node n = zero_gain_.back();
      zero_gain_.pop_back();
      if (inc_.gain(n) >= 0) return n;
      else neg_gain_.push_back(n);
    }
    Node n = neg_gain_.back();
    neg_gain_.pop_back();
    return n;
  }

  void push(Node n) {
    Weight g = inc_.gain(n);
    if (g > 0) {
      pos_gain_.push_back(n);
    }
    else if (g == 0) {
      zero_gain_.push_back(n);
    }
    else {
      neg_gain_.push_back(n);
    }
  }

  void clear() {
    pos_gain_.clear();
    zero_gain_.clear();
    neg_gain_.clear();
  }

 private:
  std::vector<Node> pos_gain_;
  std::vector<Node> zero_gain_;
  std::vector<Node> neg_gain_;
  const IncBipart &inc_;
};

void random_placement_pass(IncBipart &inc, std::minstd_rand &rgen) {
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
}

void traction_placement_pass(IncBipart &inc, std::minstd_rand &rgen) {
  // Initialize with random nodes
  std::vector<Node> nodes (inc.nodes().begin(), inc.nodes().end());
  std::shuffle(nodes.begin(), nodes.end(), rgen);

  std::size_t moves_left = 2 * inc.nNodes();

  ThresholdQueue q(inc);
  while (!inc.legal() && moves_left) {
    // Pick one node
    if (nodes.empty()) break;
    q.push(nodes.back());
    nodes.pop_back();

    // Pull other nodes with it
    // Technically possible to have quadratic complexity if the overflow moves from one side to the other
    // Hopefully very uncommon in practice
    do {
      Node n = q.pop();
      if (!inc.overflow(inc.mapping(n))) continue;
      inc.move(n, [&](Node o, Weight w) { q.push(o); });
    } while (!q.empty() && --moves_left);
  }
}

void positive_gain_pass(IncBipart &inc, std::minstd_rand &rgen) {
  std::vector<Node> nodes (inc.nodes().begin(), inc.nodes().end());
  std::shuffle(nodes.begin(), nodes.end(), rgen);

  std::vector<Node> active_set;

  for (auto n : nodes) {
    active_set.push_back(n);
  }

  while (!active_set.empty()) {
    Node n = active_set.back();
    active_set.pop_back();
    if (inc.gain(n) <= 0) continue;
    inc.tryMove(n, [&](Node o, Weight w) { active_set.push_back(o); });
  }
}

void non_negative_gain_pass(IncBipart &inc, std::minstd_rand &rgen, const int max_zero_gain_moves = 2) {
  std::vector<Node> nodes (inc.nodes().begin(), inc.nodes().end());
  std::shuffle(nodes.begin(), nodes.end(), rgen);

  ThresholdQueue q(inc);
  std::vector<int> zero_gain_moves(inc.nNodes(), 0);

  for (auto n : nodes) {
    q.push(n);
  }

  while (!q.empty()) {
    Node n = q.pop();
    Weight g = inc.gain(n);
    if (g < 0) break;
    if (g == 0 && ++zero_gain_moves[n.id] >= max_zero_gain_moves) continue;
    inc.tryMove(n, [&](Node o, int w) { q.push(o); });
  }
}

void move_all(IncBipart &inc, Edge e, const std::vector<char> &dest) {
  std::size_t i = 0;
  for (Node n : inc.nodes(e)) {
    if (inc.mapping(n) != dest[i++]) {
      inc.move(n);
    }
  }
}

void edge_centric_pass(IncBipart &inc, std::minstd_rand &rgen) {
  std::vector<Edge> edges (inc.edges().begin(), inc.edges().end());
  std::shuffle(edges.begin(), edges.end(), rgen);

  for (auto e : edges) {
    // Try to move the entire edge in both directions
    std::vector<char> initial;
    for (Node n : inc.nodes(e)) {
      initial.push_back(inc.mapping(n));
    }
    std::int64_t best_cost = inc.cost();
    int best_result = -1;
    for (int i = 0; i < 2; ++i) {
      move_all(inc, e, std::vector<char>(inc.nodes(e).size(), i));
      if (inc.cost() < best_cost) {
        best_cost = inc.cost();
        best_result = i;
      }
    }
    // Keep the best result found
    if (best_result == -1) {
      move_all(inc, e, initial);
    }
    else if (best_result == 0) {
      move_all(inc, e, std::vector<char>(inc.nodes(e).size(), 0));
    }
  }
}

void place(IncBipart &inc, std::minstd_rand &rgen) {
  traction_placement_pass(inc, rgen);
}

void optimize(IncBipart &inc, std::minstd_rand &rgen) {
  non_negative_gain_pass(inc, rgen);
}

std::vector<Mapping> place(const Problem &pb, int n_starts, std::minstd_rand &rgen) {
  std::vector<Mapping> mappings;
  for (int i = 0; i < n_starts; ++i) {
    IncBipart inc(pb);
    place(inc, rgen);
    if (!inc.legal()) continue;
    mappings.push_back(inc.mapping());
  }

  return mappings;
}

void optimize(const Problem &pb, std::vector<Mapping> &mappings, std::minstd_rand &rgen) {
  for (Mapping &m : mappings) {
    IncBipart inc(pb, m);
    optimize(inc, rgen);
    m = inc.mapping();
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
  optimize(pb, mappings, rgen);
  coarsen_recurse(pb, mappings, rgen);
}

std::vector<Mapping> solve(const Problem &pb, int n_starts) {
  std::minstd_rand rgen;

  std::vector<Mapping> mappings = place(pb, n_starts, rgen);
  optimize(pb, mappings, rgen);

  const int n_cycles = 3;
  for (int i = 0; i < n_cycles; ++i) {
    optimize_recurse(pb, mappings, rgen);
  }

  return mappings;
}

} // End namespace minipart


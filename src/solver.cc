// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#include "inc_bipart.h"
#include "coarsener.h"

#include <map>

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

class FMQueue {
 public:
  FMQueue(const IncBipart &inc)
    : inc_(inc)
    , last_bucket_(0) {
    max_gain_ = 0;
    for (Node n : inc.nodes()) {
      Weight loc_gain = 0;
      for (Edge e : inc.edges(n)) {
        loc_gain += inc.hypergraph().weight(e);
      }
      max_gain_ = std::max(loc_gain, max_gain_);
    }
    buckets_.resize(2 * max_gain_ + 1);
  }

  bool empty() {
    while (true) {
      while (last_bucket_ != 0 && buckets_[last_bucket_].empty()) {
        --last_bucket_;
      }
      if (buckets_[last_bucket_].empty()) return true;
      Node n = buckets_[last_bucket_].back();
      Weight bucket = inc_.gain(n) + max_gain_;
      if (bucket != last_bucket_) {
        buckets_[last_bucket_].pop_back();
        buckets_[bucket].push_back(n);
        last_bucket_ = std::max(last_bucket_, bucket);
      }
      else return false;
    }
  }

  Node pop() {
    empty(); // Consume all out-of-place nodes
    Node n = buckets_[last_bucket_].back();
    buckets_[last_bucket_].pop_back();
    return n;
  }

  void push(Node n) {
    Weight bucket = max_gain_ + inc_.gain(n);
    buckets_[bucket].push_back(n);
    last_bucket_ = std::max(bucket, last_bucket_);
  }

  void clear() {
    for (auto & bucket : buckets_) {
      bucket.clear();
    }
  }

 private:
  const IncBipart &inc_;
  std::vector<std::vector<Node> > buckets_;
  Weight last_bucket_;
  Weight max_gain_;
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

  // Put a limit on possible quadratic complexity: better fail early
  std::size_t moves_left = 2 * inc.nNodes();

  ThresholdQueue q(inc);
  while (!inc.legal() && moves_left) {
    // Pick one node
    if (nodes.empty()) break;
    q.push(nodes.back());
    nodes.pop_back();

    // Pull other nodes with it
    do {
      Node n = q.pop();
      if (!inc.overflow(inc.mapping(n))) continue;
      inc.move(n, [&](Node o, Weight w) { q.push(o); });
    } while (!q.empty() && --moves_left);
  }
}

template <typename Queue>
void positive_gain_pass(IncBipart &inc, std::minstd_rand &rgen) {
  std::vector<Node> nodes (inc.nodes().begin(), inc.nodes().end());
  std::shuffle(nodes.begin(), nodes.end(), rgen);

  Queue q(inc);

  for (auto n : nodes) {
    q.push_back(n);
  }

  while (!q.empty()) {
    Node n = q.pop();
    if (inc.gain(n) <= 0) continue;
    inc.tryMove(n, [&](Node o, Weight w) { q.push(o); });
  }
}

template <typename Queue>
void non_negative_gain_pass(IncBipart &inc, std::minstd_rand &rgen, const int max_zero_gain_moves = 1) {
  std::vector<Node> nodes (inc.nodes().begin(), inc.nodes().end());
  std::shuffle(nodes.begin(), nodes.end(), rgen);

  Queue q(inc);
  std::vector<int> zero_gain_moves(inc.nNodes(), 0);

  for (auto n : nodes) {
    q.push(n);
  }

  while (!q.empty()) {
    Node n = q.pop();
    Weight g = inc.gain(n);
    if (g < 0) break;
    if (g == 0 && ++zero_gain_moves[n.id] > max_zero_gain_moves) continue;
    inc.tryMove(n, [&](Node o, int w) { q.push(o); });
  }
}

template <typename Queue>
void probing_pass(IncBipart &inc, const std::vector<Node> &nodes, int max_moves_per_probe) {
  Queue q(inc);
  std::vector<Node> trail;

  for (auto n : nodes) {
    auto best_cost = inc.cost();
    int moves_left = max_moves_per_probe;

    q.push(n);
    while (!q.empty() && --moves_left) {
      Node n = q.pop();
      if (inc.tryMove(n, [&](Node o, int w) { q.push(o); })) {
        trail.push_back(n);
      }

      if (inc.cost() < best_cost) {
        trail.clear();
        moves_left = max_moves_per_probe;
        best_cost = inc.cost();
      }
    }

    for (Node n : trail) {
      inc.move(n);
    }

    q.clear();
    trail.clear();
  }
}

template <typename Queue>
void probing_pass(IncBipart &inc, std::minstd_rand &rgen, int max_moves_per_probe) {
  std::vector<Node> nodes;
  for (Node n : inc.nodes()) {
    bool frontier = false;
    for (Edge e : inc.edges(n)) {
      frontier |= inc.cut(e);
    }
    if (frontier) nodes.push_back(n);
  }
  std::shuffle(nodes.begin(), nodes.end(), rgen);
  probing_pass<Queue>(inc, nodes, max_moves_per_probe);
}

void move_all(IncBipart &inc, Edge e, const std::vector<char> &dest) {
  std::size_t i = 0;
  for (Node n : inc.nodes(e)) {
    if (inc.mapping(n) != dest[i++]) {
      inc.move(n);
    }
  }
}

void edge_centric_pass(IncBipart &inc, const std::vector<Edge> &edges) {
  for (Edge e : edges) {
    // Try to move the entire edge in both directions
    std::vector<char> initial;
    for (Node n : inc.nodes(e)) {
      initial.push_back(inc.mapping(n));
    }
    std::int64_t best_cost = inc.cost();
    int best_result = -1;
    for (int i = 0; i < 2; ++i) {
      move_all(inc, e, std::vector<char>(inc.nodes(e).size(), i));
      if (inc.legal() && inc.cost() < best_cost) {
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

void edge_centric_pass(IncBipart &inc, std::minstd_rand &rgen) {
  std::vector<Edge> edges;
  for (auto e : inc.edges()) {
    if (inc.cut(e)) edges.push_back(e);
  }
  std::shuffle(edges.begin(), edges.end(), rgen);
  edge_centric_pass(inc, edges);
}

void trySwap(IncBipart &inc, Node n1, Node n2) {
  if (n1 == n2) return;
  if (inc.mapping(n1) == inc.mapping(n2)) return;
  if (inc.gain(n1) + inc.gain(n2) <= 0) return;

  std::int64_t cost = inc.cost();
  inc.move(n1);
  inc.move(n2);
  if (inc.cost() > cost || !inc.legal()) {
    inc.move(n1);
    inc.move(n2);
  }
}

void tryAllSwaps(IncBipart &inc, const std::vector<Node> &nodes) {
  for (Node n1 : nodes) {
    if (inc.gain(n1) <= 0) continue;
    if (inc.canMove(n1)) {
      inc.move(n1);
      continue;
    }
    for (Node n2 : nodes) {
      trySwap(inc, n1, n2);
    }
  }
}

void swap_pass(IncBipart &inc, std::minstd_rand &rgen) {
  // TODO: strict runtime limit to avoid quadratic blowup

  bool positive_gain_found = false;
  for (Node n : inc.nodes()) {
    positive_gain_found |= (inc.gain(n) > 0);
  }
  if (!positive_gain_found) return;

  std::vector<Node> nodes(inc.nodes().begin(), inc.nodes().end());
  std::shuffle(nodes.begin(), nodes.end(), rgen);

  tryAllSwaps(inc, nodes);
}

void place(IncBipart &inc, std::minstd_rand &rgen) {
  random_placement_pass(inc, rgen);
}

void optimize(IncBipart &inc, std::minstd_rand &rgen) {
  non_negative_gain_pass<ThresholdQueue>(inc, rgen);
  non_negative_gain_pass<ThresholdQueue>(inc, rgen);
  probing_pass<ThresholdQueue>(inc, rgen, 5);
  //edge_centric_pass(inc, rgen);
  swap_pass(inc, rgen);
}

void place(const Problem &pb, std::vector<Mapping> &mappings, std::size_t n_starts, std::minstd_rand &rgen) {
  // Only a fixed number of times: placement may actually be difficult
  for (std::size_t i = mappings.size(); i < n_starts; ++i) {
    IncBipart inc(pb);
    place(inc, rgen);
    if (!inc.legal()) continue;
    mappings.push_back(inc.mapping());
  }
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

void sortMappings(const Problem &pb, std::vector<Mapping> &mappings) {
  std::map<std::int64_t, Mapping> sorted_mappings;
  for (Mapping &m : mappings) {
    std::int64_t cost = computeBipartCost(pb.hypergraph, m);
    sorted_mappings.emplace(cost, std::move(m));
  }
  mappings.clear();
  for (auto & p : sorted_mappings) {
    mappings.emplace_back(std::move(p.second));
  }
}

void coarsen_recurse(const Problem &pb, std::vector<Mapping> &mappings, std::minstd_rand &rgen) {
  //Coarsening coarsening = inferCoarsening(mappings);
  std::size_t target_nnodes = pb.hypergraph.nNodes() * 0.5;
  if (target_nnodes < 20) return;

  // FIXME: Keep left out mappings and insert them back?
  sortMappings(pb, mappings);
  Coarsening coarsening = selectForCoarsening(mappings, target_nnodes);

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
  mappings.clear();
  for (const Mapping &c_m : c_mappings) {
    mappings.push_back(coarsening.reverse(c_m));
  }
}

void optimize_recurse(const Problem &pb, std::vector<Mapping> &mappings, std::minstd_rand &rgen) {
  optimize(pb, mappings, rgen);
  coarsen_recurse(pb, mappings, rgen);
  optimize(pb, mappings, rgen);
}

std::vector<Mapping> solve(const Problem &pb, const SolverOptions &options) {
  std::minstd_rand rgen(options.seed == 0 ? 1 : options.seed);

  std::vector<Mapping> mappings;

  for (std::size_t i = 0; i < options.n_cycles; ++i) {
    place(pb, mappings, options.n_starts, rgen);
    optimize_recurse(pb, mappings, rgen);
  }

  return mappings;
}

} // End namespace minipart


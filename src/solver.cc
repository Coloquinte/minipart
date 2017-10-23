// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#include "coarsener.h"
#include "queues.h"

#include <map>

namespace minipart {

void legalization_pass(IncBipart &inc, std::minstd_rand &rgen) {
  if (inc.legal()) return;

  std::vector<Node> nodes (inc.nodes().begin(), inc.nodes().end());
  std::shuffle(nodes.begin(), nodes.end(), rgen);

  for (auto n : nodes) {
    if (inc.legal()) break;
    bool m = inc.mapping(n);
    if (inc.overflow(m)) inc.move(n);
  }
}

void random_placement_pass(IncBipart &inc, std::minstd_rand &rgen) {
  std::bernoulli_distribution dist;
  for (auto n : inc.nodes()) {
    if (dist(rgen)) inc.move(n);
  }

  legalization_pass(inc, rgen);
}

void traction_placement_pass(IncBipart &inc, std::minstd_rand &rgen) {
  // TODO: this is crappy: make it pull from the overflowed side

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

void greedy_pass(IncBipart &inc, std::minstd_rand &rgen, int passes=3) {
  for (int i = 0; i < passes; ++i) {
    std::vector<Node> q(inc.nodes().begin(), inc.nodes().end());
    std::shuffle(q.begin(), q.end(), rgen);
    for (Node n : q) {
      if (inc.gain(n) >= 0) {
        inc.tryMove(n);
      }
    }
  }
}

void positive_gain_pass(IncBipart &inc, std::minstd_rand &rgen) {
  std::vector<Node> q(inc.nodes().begin(), inc.nodes().end());
  std::shuffle(q.begin(), q.end(), rgen);

  while (!q.empty()) {
    Node n = q.back();
    q.pop_back();
    if (inc.gain(n) <= 0) continue;
    inc.tryMove(n, [&](Node o, Weight w) { q.push_back(o); });
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
  std::vector<char> initial;
  for (Edge e : edges) {
    // Try to move the entire edge in both directions
    for (Node n : inc.nodes(e)) {
      initial.push_back(inc.mapping(n));
    }
    std::int64_t best_cost = inc.cost();
    int best_result = -1;
    for (int i = 0; i < 2; ++i) {
      for (Node n : inc.nodes(e)) {
        if (inc.mapping(n) != i) {
          inc.move(n);
        }
      }
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
      for (Node n : inc.nodes(e)) {
        inc.move(n);
      }
    }
    initial.clear();
  }
}

void edge_centric_pass(IncBipart &inc, std::minstd_rand &rgen) {
  std::vector<Edge> edges;
  for (auto e : inc.edges()) {
    if (inc.nodes(e).size() <= 4) continue; // Mostly taken care of by other algorithms
    if (!inc.cut(e)) continue; // Only edges on the frontier
    edges.push_back(e);
  }
  std::shuffle(edges.begin(), edges.end(), rgen);
  edge_centric_pass(inc, edges);
}

void trySwap(IncBipart &inc, Node n1, Node n2) {
  if (n1 == n2) return;
  if (inc.mapping(n1) == inc.mapping(n2)) return;
  // gain(swap) <= gain(move1) + gain(move2); if non-positive, don't bother
  if (inc.gain(n1) + inc.gain(n2) <= 0) return;

  // TODO: add legality check before trying the swap
  std::int64_t cost = inc.cost();
  inc.move(n1);
  inc.move(n2);
  if (inc.cost() > cost || !inc.legal()) {
    inc.move(n1);
    inc.move(n2);
  }
}

std::vector<Node> select_biggest_nodes(const IncBipart &inc, std::size_t num) {
  std::vector<Node> ret(inc.nodes().begin(), inc.nodes().end());
  // Access the weight matrix for each resource and get the n biggest nodes
  const Matrix<Resource> &demands = inc.demands();
  assert (demands.size1() == inc.nNodes());
  // TODO: handle multiple resources
  std::sort (ret.begin(), ret.end(),
      [&](Node a, Node b) {
        return demands(a.id, 0) > demands(b.id, 0);
      });
  if (ret.size() > num) ret.resize(num);
  return ret;
}

void swap_pass(IncBipart &inc, const std::vector<Node> &nodes) {
  for (Node n1 : nodes) {
    // At least one of the nodes ought to have positive gain
    if (inc.gain(n1) <= 0) continue;
    // No point swapping if we can move
    if (inc.canMove(n1)) {
      inc.move(n1);
      continue;
    }
    for (Node n2 : nodes) {
      trySwap(inc, n1, n2);
    }
  }
}

void swap_pass(IncBipart &inc, std::minstd_rand &rgen, std::size_t num = 100) {
  std::vector<Node> nodes = select_biggest_nodes(inc, num);
  std::shuffle(nodes.begin(), nodes.end(), rgen);

  swap_pass(inc, nodes);
}

void exhaustive_pass(IncBipart &inc, const std::vector<Node> &nodes) {
  assert (nodes.size() < 64);
  std::uint64_t cnt = 1 << nodes.size();

  // Keep the flips since the last best result
  std::vector<char> history(nodes.size(), 0);
  std::int64_t best_result = inc.cost();

  std::uint64_t gray = 0;
  for (std::uint64_t i = 1; i < cnt; ++i) {
    // Convert to a Gray counter to perform only one move per iteration
    std::uint64_t cur_gray = i ^ (i >> 1);
    std::uint64_t diff = cur_gray ^ gray;
    assert ( (diff & (diff-1)) == 0); // Only one bit flipped
    gray = cur_gray;

    // Find first bit set; compiler, please optimize this i.e. count leading zeros
    std::size_t flipped = 0;
    for (std::uint64_t s = diff >> 1; s != 0; s >>= 1) {
      ++flipped;
    }

    // Make the change
    Node n = nodes[flipped];
    inc.move(n);
    history[flipped] ^= 1;

    // Reset the history
    if (inc.legal() && inc.cost() < best_result) {
      history.assign(nodes.size(), (char) 0);
      best_result = inc.cost();
    }
  }
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    if (history[i]) {
      inc.move(nodes[i]);
    }
  }
}

void exhaustive_pass(IncBipart &inc, std::minstd_rand &rgen) {
  std::vector<Node> nodes = select_biggest_nodes(inc, 10);
  std::shuffle(nodes.begin(), nodes.end(), rgen);

  exhaustive_pass(inc, nodes);
}

void place(IncBipart &inc, std::minstd_rand &rgen) {
  random_placement_pass(inc, rgen);
}

void optimize(IncBipart &inc, std::minstd_rand &rgen) {
  non_negative_gain_pass<PosQueue>(inc, rgen);
  non_negative_gain_pass<PosQueue>(inc, rgen);
  non_negative_gain_pass<PosQueue>(inc, rgen);
  non_negative_gain_pass<PosQueue>(inc, rgen);
  probing_pass<PosQueue>(inc, rgen, 5);
  edge_centric_pass(inc, rgen);
  swap_pass(inc, rgen);
  //exhaustive_pass(inc, rgen);
}

void sort_mappings(const Problem &pb, std::vector<Mapping> &mappings) {
  std::multimap<std::int64_t, Mapping> sorted_mappings;
  for (Mapping &m : mappings) {
    std::int64_t cost = computeBipartCost(pb.hypergraph, m);
    sorted_mappings.emplace(cost, std::move(m));
  }
  mappings.clear();
  for (auto & p : sorted_mappings) {
    mappings.emplace_back(std::move(p.second));
  }
}

class BipartSolver {
 public:
  BipartSolver(Problem pb, SolverOptions options, std::shared_ptr<std::vector<std::minstd_rand> > rgens, std::vector<Mapping> pool);
  BipartSolver(Problem pb, SolverOptions options);

  void init();
  void run();

  const std::vector<Mapping> &mappings() const { return solution_pool_; }

 private:
  void place();
  void optimize();
  void coarsen_recurse();

 private:
  // Problem to solve
  Problem pb_;
  // Options
  SolverOptions options_;
  // Multiple random generators for parallel algorithms
  // Shared between solvers: don't call them in parallel
  // TODO: do not share state: use Xoroshiro128+ or PCG to have
  std::shared_ptr<std::vector<std::minstd_rand> > rgens_;
  // Multiple solutions
  std::vector<Mapping> solution_pool_;
  // TODO: adaptive strategies
};

BipartSolver::BipartSolver(Problem pb, SolverOptions options)
: BipartSolver(pb, options, nullptr, std::vector<Mapping>()) {
  std::seed_seq seq = {options.seed};
  std::vector<std::uint32_t> seeds(options.n_starts);
  seq.generate(seeds.begin(), seeds.end());
  rgens_.reset(new std::vector<std::minstd_rand>());
  for (std::uint32_t s : seeds) {
    rgens_->emplace_back(s);
  }
}

BipartSolver::BipartSolver(Problem pb, SolverOptions options, std::shared_ptr<std::vector<std::minstd_rand> > rgens, std::vector<Mapping> pool)
: pb_(pb)
, options_(options)
, rgens_(rgens)
, solution_pool_(pool){
}

void BipartSolver::init() {
  place();
}

void BipartSolver::run() {
  optimize();
  coarsen_recurse();
  optimize();
}

void BipartSolver::place() {
  // Initialize the vector to access elements in parallel
  std::size_t range_begin = solution_pool_.size();
  std::size_t range_end   = options_.n_starts;
  solution_pool_.resize(range_end);

  // Reuse the same object (rather than special-purpose constructor)
  IncBipart reused(pb_);

  #pragma omp parallel for num_threads(options_.n_threads) schedule(dynamic, 1)
  for (std::size_t i = range_begin; i < range_end; ++i) {
    IncBipart inc = reused;
    minipart::place(inc, rgens_->at(i));
    // Only a fixed number of times: placement may actually be difficult and fail
    if (!inc.legal()) continue;
    solution_pool_[i] = inc.mapping();
  }
}

void BipartSolver::optimize() {
  #pragma omp parallel for num_threads(options_.n_threads) schedule(dynamic, 1)
  for (std::size_t i = 0; i < solution_pool_.size(); ++i) {
    Mapping &m = solution_pool_[i];
    IncBipart inc(pb_, m);
    minipart::optimize(inc, rgens_->at(i));
    m = inc.mapping();
  }
}

void BipartSolver::coarsen_recurse() {
  std::size_t target_nnodes = pb_.hypergraph.nNodes() * 0.5;
  if (target_nnodes < 20) return;

  // Chose a coarsening
  // FIXME: Keep left out mappings and insert them back?
  sort_mappings(pb_, solution_pool_);
  Coarsening coarsening = select_for_coarsening(solution_pool_, target_nnodes);

  // Apply the coarsening
  Problem c_pb = coarsening(pb_);
  std::vector<Mapping> c_mappings;
  for (const Mapping &m : solution_pool_) {
    auto c_m = coarsening(m);
    assert (computeBipartCost(pb_.hypergraph, m) == computeBipartCost(c_pb.hypergraph, c_m));
    c_mappings.push_back(c_m);
  }

  // Call recursively
  BipartSolver coarsened(c_pb, options_, rgens_, c_mappings);
  coarsened.run();

  // Read results back
  solution_pool_.clear();
  for (const Mapping &c_m : coarsened.mappings()) {
    auto m = coarsening.reverse(c_m);
    assert (computeBipartCost(pb_.hypergraph, m) == computeBipartCost(c_pb.hypergraph, c_m));
    solution_pool_.push_back(m);
  }
}

std::vector<Mapping> solve(const Problem &pb, const SolverOptions &options) {
  BipartSolver s(pb, options);
  for (std::size_t i = 0; i < options.n_cycles; ++i) {
    s.init();
    s.run();
  }
  auto ret = s.mappings();
  sort_mappings(pb, ret);
  return ret;
}

} // End namespace minipart


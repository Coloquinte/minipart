// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#include "queues.h"

#include <functional>

namespace minipart {

typedef std::function<void (IncBipart&, std::minstd_rand&)> Strategy;

class StrategySelector {
 public:
  StrategySelector(std::initializer_list<Strategy>, std::vector<double> weights = {});
  void run(IncBipart &, std::minstd_rand&);

 private:
  std::vector<Strategy> strategies_;
  std::vector<double> weights_;
};

StrategySelector::StrategySelector(std::initializer_list<Strategy> strategies, std::vector<double> weights)
: strategies_(strategies)
, weights_(weights) {
  weights_.resize(strategies_.size(), 1.0);
}

void StrategySelector::run(IncBipart &inc, std::minstd_rand &rgen) {
  std::discrete_distribution<> dist(weights_.begin(), weights_.end());
  int ind = dist(rgen);
  strategies_[ind](inc, rgen);
}

class StrategyComposer {
 public:
  StrategyComposer(std::initializer_list<Strategy>);
  void operator()(IncBipart &, std::minstd_rand&);

 private:
  std::vector<Strategy> strategies_;
};

StrategyComposer::StrategyComposer(std::initializer_list<Strategy> strategies)
: strategies_(strategies) {
}

void StrategyComposer::operator()(IncBipart &inc, std::minstd_rand &rgen) {
  for (auto &strategy : strategies_) {
    strategy(inc, rgen);
  }
}

bool trySwap(IncBipart &inc, Node n1, Node n2) {
  if (n1 == n2) return false;
  if (inc.mapping(n1) == inc.mapping(n2)) return false;
  // gain(swap) <= gain(move1) + gain(move2); if non-positive, don't bother
  if (inc.gain(n1) + inc.gain(n2) <= 0) return false;

  // TODO: add legality check before trying the swap
  std::int64_t cost = inc.cost();
  inc.move(n1);
  inc.move(n2);
  if (inc.cost() > cost || !inc.legal()) {
    inc.move(n1);
    inc.move(n2);
    return false;
  }
  return true;
}

void legalization_pass(IncBipart &inc, std::minstd_rand &rgen) {
  for (int i = 0; i < 4; ++i) {
    if (inc.legal()) break;

    std::vector<Node> nodes (inc.nodes().begin(), inc.nodes().end());
    std::shuffle(nodes.begin(), nodes.end(), rgen);

    for (auto n : nodes) {
      if (inc.legal()) break;
      bool m = inc.mapping(n);
      if (inc.overflow(m)) inc.move(n);
    }
  }
}

void random_placement_pass(IncBipart &inc, std::minstd_rand &rgen) {
  std::bernoulli_distribution dist;
  for (auto n : inc.nodes()) {
    if (dist(rgen)) inc.move(n);
  }

  legalization_pass(inc, rgen);
}

template <typename Queue>
void traction_placement_pass(IncBipart &inc, std::minstd_rand &rgen) {
  std::vector<Node> nodes(inc.nodes().begin(), inc.nodes().end());
  std::shuffle(nodes.begin(), nodes.end(), rgen);

  std::vector<Queue> q(2, Queue(inc));
  for (Node n : nodes) {
    q[inc.mapping(n)].push(n);
  }

  const std::size_t max_moves = 2 * inc.nNodes();

  for (std::size_t i = 0; i < max_moves; ++i) {
    if (inc.legal()) break;
    bool from = inc.overflow(true);
    while (inc.overflow(from) && !q[from].empty()) {
      Node n = q[from].pop();
      if (inc.mapping(n) == from) {
        inc.move(n, [&](Node o, int w) { q[from].push(o); });
      }
    }
  }

  legalization_pass(inc, rgen);
}

void greedy_pass(IncBipart &inc, std::minstd_rand &rgen, int passes=3) {
  for (int i = 0; i < passes; ++i) {
    std::vector<Node> q(inc.nodes().begin(), inc.nodes().end());
    std::shuffle(q.begin(), q.end(), rgen);
    for (Node n : q) {
      if (inc.gain(n) >= 0) {
        inc.tryMove(n);
      }
    ;}
  }
}

void positive_gain_pass(IncBipart &inc, std::minstd_rand &rgen) {
  std::vector<Node> q;
  for (Node n : inc.nodes()) {
    if (inc.gain(n) > 0) q.push_back(n);
  }
  std::shuffle(q.begin(), q.end(), rgen);

  while (!q.empty()) {
    Node n = q.back();
    q.pop_back();
    if (inc.gain(n) <= 0) continue;
    inc.tryMove(n, [&](Node o, Weight w) { q.push_back(o); });
  }
}

template <typename Queue>
void non_negative_gain_pass(IncBipart &inc, std::minstd_rand &rgen) {
  const int max_zero_gain_moves = 1;
  std::vector<Node> nodes;
  for (Node n : inc.nodes()) {
    if (inc.gain(n) >= 0) nodes.push_back(n);
  }
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
void dual_queue_pass(IncBipart &inc, std::minstd_rand &rgen) {
  const int max_zero_gain_moves = 1;
  std::vector<Node> nodes;
  for (Node n : inc.nodes()) {
    if (inc.gain(n) >= 0) nodes.push_back(n);
  }
  std::shuffle(nodes.begin(), nodes.end(), rgen);

  std::vector<Queue> q(2, Queue(inc));
  for (Node n : nodes) {
    q[inc.mapping(n)].push(n);
  }
  std::vector<int> zero_gain_moves(inc.nNodes(), 0);

  bool current = false;
  while (!q[0].empty() || !q[1].empty()) {
    while (!q[current].empty()) {
      Node n = q[current].pop();
      Weight g = inc.gain(n);
      if (g < 0) continue;
      if (g == 0 && ++zero_gain_moves[n.id] > max_zero_gain_moves) continue;
      if (inc.mapping(n) == current && inc.canMove(n)) {
        inc.move(n, [&](Node o, int w) { q[current].push(o); });
      }
    }
    current = !current;
  }
}

template <typename Queue>
class HybridPassFunctor {
 public:
  HybridPassFunctor(IncBipart &inc, std::minstd_rand &rgen, int max_moves_per_probe);
  void run();

 private:
  bool empty();
  bool attempt_move();
  bool attempt_swap();
  void drop();

 private:
  IncBipart &inc_;
  std::minstd_rand &rgen_;
  std::vector<int> zero_gain_moves_;
  int max_zero_gain_moves_;
  std::vector<Queue> q_;
  bool current_;
};

template <typename Queue>
HybridPassFunctor<Queue>::HybridPassFunctor(IncBipart &inc, std::minstd_rand &rgen, int max_zero_gain_moves)
: inc_(inc)
, rgen_(rgen)
, max_zero_gain_moves_(max_zero_gain_moves) {
  zero_gain_moves_.assign(inc.nNodes(), 0);
  q_ = std::vector<Queue> (2, Queue(inc));
  current_ = false;

  std::vector<Node> nodes;
  for (Node n : inc.nodes()) {
    if (inc.gain(n) >= 0) nodes.push_back(n);
  }
  std::shuffle(nodes.begin(), nodes.end(), rgen);

  for (Node n : nodes) {
    q_[inc.mapping(n)].push(n);
  }
}

template <typename Queue>
bool HybridPassFunctor<Queue>::empty() {
  return q_[0].empty() && q_[1].empty();
}

template <typename Queue>
void HybridPassFunctor<Queue>::run() {
  while (!empty()) {
    bool move_success = attempt_move();
    if (move_success) continue;

    bool swap_success = attempt_swap();
    if (swap_success) continue;

    drop();
  }
}

template <typename Queue>
bool HybridPassFunctor<Queue>::attempt_move() {
  bool success = false;
  // Consume all moves on both queues
  Queue &q = q_[current_];
  while (!q.empty()) {
    Node n = q.pop();
    Weight g = inc_.gain(n);
    if (g < 0) continue;
    if (g == 0 && ++zero_gain_moves_[n.id] > max_zero_gain_moves_) continue;
    if (inc_.mapping(n) == current_ && inc_.canMove(n)) {
      inc_.move(n, [&](Node o, int w) { q.push(o); });
      success = true;
    }
    else {
      // Keep it for later
      q.push(n);
      break;
    }
  }
  current_ = !current_;
  return success;
}

template <typename Queue>
bool HybridPassFunctor<Queue>::attempt_swap() {
  if (!q_[0].empty() && !q_[1].empty()) {
    Node n0 = q_[0].pop();
    Node n1 = q_[1].pop();
    bool success = trySwap(inc_, n0, n1);
    if (!success) {
      q_[0].push(n0);
      q_[1].push(n1);
    }
  }
  return false;
}

template <typename Queue>
void HybridPassFunctor<Queue>::drop() {
  // Drop an unmovable node from the queue
  bool empty0 = q_[0].empty();
  bool empty1 = q_[1].empty();
  if (!empty0 && !empty1) {
    std::bernoulli_distribution dist;
    bool side = dist(rgen_);
    q_[side].pop();
  }
  else if (!empty0) q_[0].pop();
  else if (!empty1) q_[1].pop();
}

template <typename Queue>
void hybrid_pass(IncBipart &inc, std::minstd_rand &rgen) {
  const int max_zero_gain_moves = 1;
  HybridPassFunctor<Queue> pass(inc, rgen, max_zero_gain_moves);
  pass.run();
}

template <typename Queue>
void probing_pass_helper(IncBipart &inc, const std::vector<Node> &nodes, int max_moves_per_probe) {
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
void probing_pass(IncBipart &inc, std::minstd_rand &rgen) {
  int max_moves_per_probe = 5;
  std::vector<Node> nodes;
  for (Node n : inc.nodes()) {
    bool frontier = false;
    for (Edge e : inc.edges(n)) {
      frontier |= inc.cut(e);
    }
    if (frontier) nodes.push_back(n);
  }
  std::shuffle(nodes.begin(), nodes.end(), rgen);
  probing_pass_helper<Queue>(inc, nodes, max_moves_per_probe);
}

void move_all(IncBipart &inc, Edge e, const std::vector<char> &dest) {
  std::size_t i = 0;
  for (Node n : inc.nodes(e)) {
    if (inc.mapping(n) != dest[i++]) {
      inc.move(n);
    }
  }
}

void edge_centric_pass_helper(IncBipart &inc, const std::vector<Edge> &edges) {
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
  edge_centric_pass_helper(inc, edges);
}

std::vector<Node> select_biggest_nodes(const IncBipart &inc, std::size_t num) {
  // Access the weight matrix for each resource and get the n biggest nodes
  typedef std::pair<Node, Resource> P;
  const Matrix<Resource> &demands = inc.demands();

  // TODO: handle multiple resources
  std::vector<P> nodes;
  for (Node n : inc.nodes()) {
    nodes.emplace_back(n, demands(n.id, 0));
  }

  std::sort (nodes.begin(), nodes.end(),
      [](P a, P b) { return a.second > b.second; });
  if (nodes.size() > num) nodes.resize(num);

  std::vector<Node> ret;
  ret.reserve(num);
  for (P p : nodes) {
    ret.push_back(p.first);
  }
  return ret;
}

void swap_pass_helper(IncBipart &inc, const std::vector<Node> &nodes) {
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

void swap_pass(IncBipart &inc, std::minstd_rand &rgen) {
  std::size_t num = 100;
  std::vector<Node> nodes = select_biggest_nodes(inc, num);
  std::shuffle(nodes.begin(), nodes.end(), rgen);

  swap_pass_helper(inc, nodes);
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

void empty_pass(IncBipart &inc, std::minstd_rand &rgen) {
}

void place(IncBipart &inc, std::minstd_rand &rgen, SolverOptions options) {
  StrategySelector sel ({
      random_placement_pass
    , traction_placement_pass<ThresholdQueue<false> >
    , traction_placement_pass<ThresholdQueue<true> >
    , traction_placement_pass<BasicQueue<false> >
    , traction_placement_pass<BasicQueue<true> >
    , traction_placement_pass<FMQueue<false> >
  }, options.place_strategies);
  sel.run(inc, rgen);
}

void optimize(IncBipart &inc, std::minstd_rand &rgen, SolverOptions options) {
  StrategySelector sel_descent ({
      hybrid_pass<PosQueue<true> >
    , hybrid_pass<PosQueue<false> >
    , hybrid_pass<FMQueue<true> >
    , hybrid_pass<FMQueue<false> >
    , non_negative_gain_pass<PosQueue<false> >
    , non_negative_gain_pass<FMQueue<false> >
  });

  StrategySelector sel_special ({
      StrategyComposer({
        edge_centric_pass
      , swap_pass})
    , StrategyComposer({
        swap_pass
      , edge_centric_pass})
    , probing_pass<PosQueue<> >
    , probing_pass<ThresholdQueue<> >
    , empty_pass
  });

  sel_descent.run(inc, rgen);
  sel_special.run(inc, rgen);
  sel_descent.run(inc, rgen);
}

} // End namespace minipart


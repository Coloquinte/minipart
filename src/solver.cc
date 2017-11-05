// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#include "coarsener.h"
#include "inc_bipart.h"

#include <iostream>

namespace minipart {

class BipartSolver {
 public:
  BipartSolver(Problem pb, SolverOptions options, std::shared_ptr<std::vector<std::minstd_rand> > rgens, std::vector<Mapping> pool, std::size_t level);
  BipartSolver(Problem pb, SolverOptions options);

  void place();
  void run();
  const Mapping &solution();

  const std::vector<Mapping> &mappings() const { return solution_pool_; }

 private:
  void optimize();
  void coarsen_recurse();
  void report();
  void sort_pool();

 private:
  // Problem to solve
  Problem pb_;
  // Options
  SolverOptions options_;
  // Coarsening level
  std::size_t level_;
  // Multiple random generators for parallel algorithms
  // Shared between solvers: don't call them in parallel
  // TODO: do not share state: use Xoroshiro128+ or PCG to have
  std::shared_ptr<std::vector<std::minstd_rand> > rgens_;
  // Multiple solutions
  std::vector<Mapping> solution_pool_;
  // TODO: adaptive strategies
};

BipartSolver::BipartSolver(Problem pb, SolverOptions options)
: BipartSolver(pb, options, nullptr, std::vector<Mapping>(), 0) {
  std::seed_seq seq = {options.seed};
  std::vector<std::uint32_t> seeds(options.n_starts);
  seq.generate(seeds.begin(), seeds.end());
  rgens_.reset(new std::vector<std::minstd_rand>());
  for (std::uint32_t s : seeds) {
    rgens_->emplace_back(s);
  }
}

BipartSolver::BipartSolver(Problem pb, SolverOptions options, std::shared_ptr<std::vector<std::minstd_rand> > rgens, std::vector<Mapping> pool, std::size_t level)
: pb_(pb)
, options_(options)
, level_(level)
, rgens_(rgens)
, solution_pool_(pool) {
}

void BipartSolver::run() {
  optimize();
  report();
  coarsen_recurse();
  optimize();
  report();
}

void BipartSolver::report() {
  if (options_.verbosity <= 1) return;

  for (std::size_t i = 0; i < level_; ++i) std::cout << "  ";
  std::cout << pb_.hypergraph.nNodes() << " nodes, ";
  std::cout << solution_pool_.size() << " solutions, ";
  std::int64_t min_cost = std::numeric_limits<std::int64_t>::max();
  double avg_cost = 0.0;
  for (const Mapping &m : solution_pool_) {
    std::int64_t cost = computeCostBipart(pb_.hypergraph, m);
    min_cost = std::min(min_cost, cost);
    avg_cost += cost;
  }
  avg_cost /= solution_pool_.size();
  std::cout << avg_cost << " average cost, ";
  std::cout << min_cost << " minimum cost" << std::endl;
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
    minipart::place(inc, rgens_->at(i), options_);
    // Only a fixed number of times: placement may actually be difficult and fail
    if (inc.legal()) {
      solution_pool_[i] = inc.mapping();
    }
  }

  // Erase positions that haven't been filled
  solution_pool_.erase(
    std::remove_if(solution_pool_.begin(), solution_pool_.end()
      , [](const Mapping &m) { return m.nNodes() == 0; }),
      solution_pool_.end()
  );
}

void BipartSolver::optimize() {
  #pragma omp parallel for num_threads(options_.n_threads) schedule(dynamic, 1)
  for (std::size_t i = 0; i < solution_pool_.size(); ++i) {
    Mapping &m = solution_pool_[i];
    IncBipart inc(pb_, m);
    minipart::optimize(inc, rgens_->at(i), options_);
    m = inc.mapping();
  }
}

void BipartSolver::coarsen_recurse() {
  std::size_t target_n_nodes = pb_.hypergraph.nNodes() * 0.8;
  if (target_n_nodes < 20) return;

  sort_pool();
  auto next_pools = select_pool_coarsenings(pb_, solution_pool_, target_n_nodes, rgens_->at(0));
  solution_pool_.clear();

  for (std::pair<Coarsening, std::vector<Mapping> > &coarsening_pool : next_pools) {
    const Coarsening &coarsening = coarsening_pool.first;
    std::vector<Mapping> &selected = coarsening_pool.second;
    // Apply the coarsening
    Problem c_pb = coarsening(pb_);
    std::vector<Mapping> c_mappings;
    for (Mapping &m : selected) {
      // TODO: do those in parallel
      c_mappings.emplace_back(coarsening(m));
      assert (computeCostBipart(pb_.hypergraph, m) == computeCostBipart(c_pb.hypergraph, c_mappings.back()));
      m = Mapping(); // Release memory early
    }

    // Call recursively
    BipartSolver coarsened(c_pb, options_, rgens_, c_mappings, level_+1);
    coarsened.run();

    // Read results back
    for (const Mapping &c_m : coarsened.solution_pool_) {
      solution_pool_.emplace_back(coarsening.reverse(c_m));
      assert (computeCostBipart(pb_.hypergraph, solution_pool_.back()) == computeCostBipart(c_pb.hypergraph, c_m));
    }
  }
}

void BipartSolver::sort_pool() {
  typedef std::pair<std::int64_t, Mapping> CM;
  std::vector<CM> cost_to_mapping;
  for (const Mapping & m : solution_pool_) {
    cost_to_mapping.emplace_back(computeCostBipart(pb_.hypergraph, m), m);
  }
  std::sort (cost_to_mapping.begin(), cost_to_mapping.end(), [](const CM &a, const CM &b) { return a.first < b.first; });
  solution_pool_.clear();
  for (CM &c : cost_to_mapping) {
    solution_pool_.emplace_back(c.second);
  }
}

const Mapping &BipartSolver::solution() {
  if (solution_pool_.empty()) throw std::runtime_error("No solution found");
  sort_pool();
  return solution_pool_.front();
}

Mapping solve(const Problem &pb, const SolverOptions &options) {
  BipartSolver s(pb, options);
  s.place();
  for (std::size_t i = 0; i < options.n_cycles; ++i) {
    s.run();
  }
  return s.solution();
}

} // End namespace minipart


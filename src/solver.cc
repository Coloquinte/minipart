// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#include "coarsener.h"
#include "inc_bipart.h"

#include <map>

namespace minipart {

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


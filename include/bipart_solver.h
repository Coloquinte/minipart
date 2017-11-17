// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#pragma once

#include "solver.h"

namespace minipart {

class BipartSolver {
 public:
  BipartSolver(Problem pb, SolverOptions options);

  void run();

  Mapping solution();
  bool success() const;

  const std::vector<Mapping> &mappings() const { return solution_pool_; }

 private:
  BipartSolver(Problem pb, SolverOptions options, std::shared_ptr<std::vector<std::minstd_rand> > rgens, std::vector<Mapping> pool, std::size_t level);
  void place();
  void run_optim();
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

} // End namespace minipart

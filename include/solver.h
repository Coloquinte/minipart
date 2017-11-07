// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#pragma once

#include "hypergraph.h"

namespace minipart {

struct Problem {
  Hypergraph hypergraph;
  Matrix<Resource> demands;
  Matrix<Resource> capacities;

  void checkConsistency() const;
};

struct SolverOptions {
  std::size_t n_starts;
  std::size_t n_cycles;
  std::size_t n_threads;
  std::size_t seed;
  std::size_t verbosity;

  std::vector<double> place_strategies;
  std::vector<double> search_strategies;
};

// Main solving function
Mapping solve(const Problem &pb, const SolverOptions &options);
Mapping bipart_solve(const Problem &pb, const SolverOptions &options);

}  // End namespace minipart


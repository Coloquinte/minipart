// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#pragma once

#include "hypergraph.h"

namespace minipart {

struct Problem {
  Hypergraph hypergraph;
  Matrix<Resource> demands;
  Matrix<Resource> capacities;
};

struct SolverOptions {
  std::size_t n_starts;
  std::size_t n_cycles;
  std::size_t n_threads;
  std::size_t seed;
  std::size_t verbosity;
};

// Main solving function
Mapping solve(const Problem &pb, const SolverOptions &options);

}  // End namespace minipart


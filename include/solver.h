// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#pragma once

#include "hypergraph.h"

namespace minipart {

struct Problem {
  Hypergraph hypergraph;
  Matrix<Resource> demands;
  Matrix<Resource> capacities;
};

std::int64_t computeBipartCost(const Hypergraph&, const Mapping&);

struct SolverOptions {
  int n_starts;
  std::size_t seed;
};

// Main solving function
std::vector<Mapping> solve(const Problem &pb, const SolverOptions &options);

}  // End namespace minipart


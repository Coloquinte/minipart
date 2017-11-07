// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#include "solver.h"

namespace minipart {

Mapping solve(const Problem &pb, const SolverOptions &options) {
  // Recursive bipartitioning
  return bipart_solve(pb, options);
}

} // End namespace minipart


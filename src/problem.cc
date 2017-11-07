// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#include "solver.h"

namespace minipart {

void Problem::checkConsistency() const {
  if (demands.size1() != hypergraph.nNodes()) abort();
  if (demands.size2() != capacities.size2()) abort();
}

}  // End namespace minipart


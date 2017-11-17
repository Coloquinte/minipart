// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#include "solver.h"

namespace minipart {

void Problem::check_consistency() const {
  if (demands.size1() != hypergraph.nNodes()) abort();
  if (demands.size2() != capacities.size2()) abort();
}

bool Problem::is_legal(const Mapping &m) const {
  if (m.nNodes() != nNodes()) return false;

  Matrix<Resource> usage = boost::numeric::ublas::zero_matrix<Resource>(nParts(), nResources());

  for (Node n : nodes()) {
    for (std::size_t j = 0; j < nResources(); ++j) {
      if (m[n].id >= nParts()) return false;
      usage(m[n].id, j) += demands(n.id, j);
    }
  }

  for (std::size_t i = 0; i < nParts(); ++i) {
    for (std::size_t j = 0; j < nResources(); ++j) {
      if (usage(i, j) > capacities(i, j)) {
        return false;
      }
    }
  }

  return true;
}

}  // End namespace minipart


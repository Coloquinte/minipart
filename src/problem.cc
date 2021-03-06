// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#include "solver.h"

namespace minipart {

namespace {
std::vector<Resource> computeTotal(const Matrix<Resource> &mat) {
  std::vector<int> totals(mat.size2(), 0);
  for (std::size_t i = 0; i < mat.size1(); ++i) {
    for (std::size_t j = 0; j < mat.size2(); ++j) {
      totals[j] += mat(i, j);
    }
  }
  return totals;
}
}

std::vector<Resource> Problem::getTotalDemands() const {
  return computeTotal(demands);
}
std::vector<Resource> Problem::getTotalCapacities() const {
  return computeTotal(capacities);
}

void Problem::check_consistency() const {
  if (demands.size1() != hypergraph.nNodes())  throw std::runtime_error("Inconsistent number of nodes");
  if (demands.size2() != capacities.size2()) throw std::runtime_error("Inconsistent number of resources");
  hypergraph.checkConsistency();
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


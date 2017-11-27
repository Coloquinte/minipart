// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#pragma once

#include "hypergraph.h"

namespace minipart {

struct Problem {
  Hypergraph hypergraph;
  Matrix<Resource> demands;
  Matrix<Resource> capacities;

  std::size_t nNodes() const { return hypergraph.nNodes(); }
  std::size_t nParts() const { return capacities.size1(); }
  std::size_t nResources() const { return capacities.size2(); }

  Range<Node> nodes() const { return hypergraph.nodes(); }
  Range<Part> parts() const { return Range<Part>(Part(0), Part(nParts())); }

  std::vector<Resource> getTotalDemands() const;
  std::vector<Resource> getTotalCapacities() const;

  void check_consistency() const;

  bool is_legal(const Mapping &) const;
};

}  // End namespace minipart



// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#pragma once

#include "solver.h"

namespace minipart {

class Coarsening {
 public:
  Coarsening(std::size_t n_in=0, std::size_t n_out=1) : m_(n_in, Node(0)), n_out_(n_out) {}

  std::size_t nNodesIn() const { return m_.size(); }
  std::size_t nNodesOut() const { return n_out_; }

  Node  operator[](Node n) const { return m_[n.id]; }
  Node& operator[](Node n)       { return m_[n.id]; }

  // Apply a coarsening
  Mapping operator() (const Mapping &m) const;
  Problem operator() (const Problem &p) const;

  // Reverse a coarsening
  Mapping reverse (const Mapping &m) const;

  void checkConsistency() const;

 private:
  std::vector<Node> m_;
  std::size_t n_out_;
};

// Create new pool levels from an existing pool
std::vector<std::pair<Coarsening, std::vector<Mapping> > > select_pool_coarsenings(const Problem &pb, const std::vector<Mapping> &pool, std::size_t target_n_nodes, std::minstd_rand &rgen);

}  // End namespace minipart


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

Coarsening inferCoarsening(const std::vector<Mapping> &mappings);
Coarsening selectForCoarsening(std::vector<Mapping> &mappings, std::size_t target_nnodes);

}  // End namespace minipart


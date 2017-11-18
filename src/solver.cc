// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#include "solver.h"
#include "bipart_solver.h"

#include <unordered_map>
#include <bitset>
#include <iostream>

namespace minipart {

class BisectionState {
 public:
  BisectionState(const Problem &pb, const SolverOptions &options);
  void run();
  Mapping mapping() const;

  void checkConsistency() const;

 private:
  struct BisectionProblem {
    Problem problem;
    std::vector<Part> parts_lft;
    std::vector<Part> parts_rgt;
  };
  BisectionProblem getBisectionProblem(Index i) const;
  void bisect();
  Index nCurrentParts() const { return subparts_.size(); }

 private:
  const Problem &pb_;
  SolverOptions options_;
  std::vector<std::vector<Part> > subparts_;
  std::vector<std::vector<Node> > part_to_nodes_;
};

BisectionState::BisectionState(const Problem &pb, const SolverOptions &options)
: pb_(pb)
, options_(options) {
  subparts_.emplace_back(pb.parts().begin(), pb.parts().end());
  part_to_nodes_.emplace_back(pb.nodes().begin(), pb.nodes().end());
}

inline bool isEdgeCut(const Hypergraph &h, const Mapping &m, Edge e) {
  // TODO: implement large version
  std::bitset<64> used;
  for (auto n : h.nodes(e)) {
    used.set(m[n].id);
  }
  return used.count() >= 2;
}

void BisectionState::run() {
  while (nCurrentParts() < pb_.nParts()) {
    if (options_.verbosity >= 2) {
      std::cout << "Bisection from " << nCurrentParts() << " partitions" << std::endl;
    }
    bisect();
    checkConsistency();
  }
}

void BisectionState::bisect() {
  std::vector<std::vector<Part> > next_subparts;
  std::vector<std::vector<Node> > next_part_to_nodes;

  for (Index i = 0; i < nCurrentParts(); ++i) {
    if (subparts_[i].size() == 1) {
      // Trivial case
      next_subparts.emplace_back(subparts_[i]);
      next_part_to_nodes.emplace_back(part_to_nodes_[i]);
      continue;
    }

    // Extract the problem itself
    BisectionProblem bisection = getBisectionProblem(i);

    // Run bipartitioning
    BipartSolver s(bisection.problem, options_);
    s.run();
    if (!s.success()) {
      std::cerr << "No bipartitioning solution found" << std::endl;
      exit(1);
    }
    Mapping sol = s.solution();

    std::vector<Node> nodes_lft, nodes_rgt;
    for (Index j = 0; j < part_to_nodes_[i].size(); ++j) {
      Node n = part_to_nodes_[i][j];
      Node o = Node(j);
      Part p = sol[o];
      if (p.id) nodes_rgt.push_back(n);
      else      nodes_lft.push_back(n);
    }
    next_part_to_nodes.push_back(nodes_lft);
    next_part_to_nodes.push_back(nodes_rgt);
    next_subparts.push_back(bisection.parts_lft);
    next_subparts.push_back(bisection.parts_rgt);
  }

  subparts_ = next_subparts;
  part_to_nodes_ = next_part_to_nodes;
}

BisectionState::BisectionProblem BisectionState::getBisectionProblem(Index i) const {
  BisectionProblem ret;
  assert (subparts_[i].size() >= 2);

  // Setup the subparts
  auto middle = subparts_[i].begin() + (subparts_[i].size() + 1) / 2;
  ret.parts_lft.assign(subparts_[i].begin(), middle);
  ret.parts_rgt.assign(middle, subparts_[i].end());
  assert (!ret.parts_lft.empty() && !ret.parts_rgt.empty());

  // Setup the demands
  ret.problem.demands = Matrix<Resource>(part_to_nodes_[i].size(), pb_.nResources());
  for (Index j = 0; j < part_to_nodes_[i].size(); ++j) {
    Node n = part_to_nodes_[i][j];
    for (Index k = 0; k < pb_.nResources(); ++k) {
      ret.problem.demands(j, k) = pb_.demands(n.id, k);
    }
  }

  // Setup the capacities
  ret.problem.capacities = boost::numeric::ublas::zero_matrix<Resource>(2, pb_.nResources());
  for (Part p : ret.parts_lft) {
    for (Index j = 0; j < pb_.nResources(); ++j) {
      ret.problem.capacities(0, j) += pb_.capacities(p.id, j);
    }
  }
  for (Part p : ret.parts_rgt) {
    for (Index j = 0; j < pb_.nResources(); ++j) {
      ret.problem.capacities(1, j) += pb_.capacities(p.id, j);
    }
  }

  // Setup the hypergraph
  std::unordered_map<Index, std::vector<Node> > edge_to_nodes;
  for (Index j = 0; j < part_to_nodes_[i].size(); ++j) {
    Node n = part_to_nodes_[i][j];
    for (Edge e : pb_.hypergraph.edges(n)) {
      edge_to_nodes[e.id].push_back(Node(j));
    }
  }

  HypergraphBuilder b(part_to_nodes_[i].size());
  for (auto const &p : edge_to_nodes) {
    Edge e(p.first);
    if (p.second.size() == 1) continue;

    // If the cost function is just the cut, only use edges that are not cut
    if (!options_.soed_objective && p.second.size() != pb_.hypergraph.nodes(e).size()) continue;

    b.addEdge(p.second.begin(), p.second.end(), pb_.hypergraph.weight(e));
  }
  ret.problem.hypergraph = b;

  return ret;
}

void BisectionState::checkConsistency() const {
  // TODO
}

Mapping BisectionState::mapping() const {
  Mapping m(pb_.nNodes());

  for (Index i = 0; i < nCurrentParts(); ++i) {
    for (Node n : part_to_nodes_[i]) {
      m[n] = Part(i);
    }
  }

  return m;
}

Mapping solve(const Problem &pb, const SolverOptions &options) {
  BisectionState state(pb, options);
  state.run();

  return state.mapping();
}

} // End namespace minipart


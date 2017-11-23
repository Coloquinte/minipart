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

  Index nCurrentParts() const { return subparts_.size(); }
  void checkConsistency() const;

 private:
  // One round of bisection
  void bisect();
  // Just setup the state
  void doBisection();
  // Obtain a legal initial solution
  void legalizeBisection();
  // Run the real thing
  void optimizeBisection();
  // Redo partitioning
  void refineBisection();

  struct BisectionProblem {
    Problem problem;
    std::vector<Node> nodes;
    Mapping mapping;
  };
  BisectionProblem getBisectionProblem(Index i, Index j) const;
  bool legalizePartitions(const std::vector<Index> &vecs);
  void redoBisection(Index i, Index j);

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
  doBisection();
  legalizeBisection();
  optimizeBisection();
  refineBisection();
}

void BisectionState::doBisection() {
  std::vector<std::vector<Part> > next_subparts;
  std::vector<std::vector<Node> > next_part_to_nodes;

  // Just setup the bisection (illegal solution)
  for (Index i = 0; i < nCurrentParts(); ++i) {
    if (subparts_[i].size() == 1) continue;
    auto middle = subparts_[i].begin() + (subparts_[i].size() + 1) / 2;
    next_subparts.emplace_back(subparts_[i].begin(), middle);
    next_subparts.emplace_back(middle, subparts_[i].end());
    next_part_to_nodes.emplace_back(part_to_nodes_[i]);
    next_part_to_nodes.emplace_back();
  }
  for (Index i = 0; i < nCurrentParts(); ++i) {
    if (subparts_[i].size() != 1) continue;
    next_subparts.emplace_back(subparts_[i]);
    next_part_to_nodes.emplace_back(part_to_nodes_[i]);
  }

  subparts_ = next_subparts;
  part_to_nodes_ = next_part_to_nodes;
}

void BisectionState::legalizeBisection() {
  // Run a bin-packing pass on all pairs to get a legal placement
  // TODO
}

void BisectionState::optimizeBisection() {
  for (Index i = 0; i+1 < nCurrentParts(); i += 2) {
    redoBisection(i, i+1);
  }
}

void BisectionState::refineBisection() {
  // TODO
}

void BisectionState::redoBisection(Index i, Index j) {
  // Extract the problem itself
  BisectionProblem bisection = getBisectionProblem(i, j);

  // Run bipartitioning
  BipartSolver s(bisection.problem, options_);
  // TODO: use initial legal solution
  //s.addInitialSolution(bisection.mapping);
  s.run();
  if (!s.success()) {
    std::cerr << "No bipartitioning solution found" << std::endl;
    exit(1);
  }
  Mapping sol = s.solution();

  // Clear the current state
  part_to_nodes_[i].clear();
  part_to_nodes_[j].clear();

  // Write a new state
  for (Index k = 0; k < bisection.nodes.size(); ++k) {
    Node n = bisection.nodes[k];
    Node o = Node(k);
    Part p = sol[o];
    if (!p.id) part_to_nodes_[i].push_back(n);
    else       part_to_nodes_[j].push_back(n);
  }
}

bool BisectionState::legalizePartitions(const std::vector<Index> &parts) {
  // Get the nodes
  // Get the resources
  // Run a greedy bin packing, biggest nodes first
  // Return success
  return true;
}

BisectionState::BisectionProblem BisectionState::getBisectionProblem(Index i, Index j) const {
  BisectionProblem ret;
  assert (!subparts_[i].empty() && !subparts_[j].empty());

  // Setup the nodes
  for (Node n : part_to_nodes_[i]) ret.nodes.push_back(n);
  for (Node n : part_to_nodes_[j]) ret.nodes.push_back(n);

  // Setup the mapping
  ret.mapping = Mapping(ret.nodes.size());
  Index m_index = 0;
  for (Index k = 0; k < part_to_nodes_[i].size(); ++k) ret.mapping[Node(m_index++)] = Part(0);
  for (Index k = 0; k < part_to_nodes_[j].size(); ++k) ret.mapping[Node(m_index++)] = Part(1);

  // Setup the demands
  ret.problem.demands = Matrix<Resource>(ret.nodes.size(), pb_.nResources());
  for (Index k = 0; k < ret.nodes.size(); ++k) {
    Node n = ret.nodes[k];
    for (Index l = 0; l < pb_.nResources(); ++l) {
      ret.problem.demands(k, l) = pb_.demands(n.id, l);
    }
  }

  // Setup the capacities
  ret.problem.capacities = boost::numeric::ublas::zero_matrix<Resource>(2, pb_.nResources());
  for (Part p : subparts_[i]) {
    for (Index j = 0; j < pb_.nResources(); ++j) {
      ret.problem.capacities(0, j) += pb_.capacities(p.id, j);
    }
  }
  for (Part p : subparts_[j]) {
    for (Index j = 0; j < pb_.nResources(); ++j) {
      ret.problem.capacities(1, j) += pb_.capacities(p.id, j);
    }
  }

  // TODO: Decrease margins for problems with multiple partitions

  // Setup the hypergraph
  std::unordered_map<Index, std::vector<Node> > edge_to_nodes;
  Index ind = 0;
  for (Index k = 0; k < ret.nodes.size(); ++k) {
    Node n = ret.nodes[k];
    Node o(ind++);
    for (Edge e : pb_.hypergraph.edges(n)) {
      edge_to_nodes[e.id].push_back(o);
    }
  }

  HypergraphBuilder b(ret.nodes.size());
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


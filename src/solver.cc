// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#include "solver.h"
#include "bipart_solver.h"

#include <bitset>
#include <iostream>

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

class BisectionState {
 public:
  BisectionState(const Problem &pb);
  void run(const SolverOptions&);
  const Mapping &mapping() const { return m_; }

  void bisect(const SolverOptions&);
  Index nCurrentParts() const { return subparts_.size(); }

  void checkConsistency() const;

 private:
  std::vector<Problem> extractBisectionProblems(bool degree=false) const;

 private:
  const Problem &pb_;
  Mapping m_;
  std::vector<std::vector<Part> > subparts_;
};

BisectionState::BisectionState(const Problem &pb) 
: pb_(pb)
, m_(pb.nNodes()) {
  std::vector<Part> parts;
  for (std::size_t i = 0; i < pb.nParts(); ++i) {
    parts.emplace_back(i);
  }
  subparts_.push_back(parts);
}

inline bool isEdgeCut(const Hypergraph &h, const Mapping &m, Edge e) {
  // TODO: implement large version
  std::bitset<64> used;
  for (auto n : h.nodes(e)) {
    used.set(m[n].id);
  }
  return used.count() >= 2;
}

void BisectionState::run(const SolverOptions &options) {
  while (nCurrentParts() < pb_.nParts()) {
    if (options.verbosity >= 2) {
      std::cout << "Bisection from " << nCurrentParts() << " partitions" << std::endl;
    }
    bisect(options);
    checkConsistency();
  }
}

void BisectionState::bisect(const SolverOptions &options) {
  // TODO: refactor this monstruosity
  // One problem for each subpart
  std::vector<Problem> problems(subparts_.size());
  const Hypergraph &h = pb_.hypergraph;

  std::vector<std::vector<Node> > part_to_nodes(subparts_.size());
  std::vector<Node> node_mapping(pb_.nNodes());
  for (Node n : pb_.nodes()) {
    Part p = m_[n];
    node_mapping[n.id] = Node(part_to_nodes[p.id].size());
    part_to_nodes[p.id].push_back(n);
  }

  // Create the hypergraph
  std::vector<HypergraphBuilder> hypergraphs;
  for (const auto &nodes : part_to_nodes) {
    hypergraphs.emplace_back(nodes.size());
  }

  for (Edge e : h.edges()) {
    // Cut cost --> export edges that are completely in one subpart
    // TODO: Edge degree cost --> export edges to all subparts where they have multiple pins
    if (isEdgeCut(h, m_, e)) continue;
    if (h.nodes(e).empty()) continue;
    Part p = m_[h.nodes(e)[0]];
    std::vector<Node> edge_nodes(h.nodes(e).begin(), h.nodes(e).end());
    for (Node &n : edge_nodes) {
      n = node_mapping[n.id];
    }
    hypergraphs[p.id].addEdge(edge_nodes.begin(), edge_nodes.end(), h.weight(e));
  }

  for (Index i = 0; i < nCurrentParts(); ++i) {
    if (subparts_[i].size() == 1) continue;
    // Reserve size
    problems[i].hypergraph = hypergraphs[i];
    problems[i].demands = Matrix<Resource>(part_to_nodes[i].size(), pb_.nResources());
    problems[i].capacities = boost::numeric::ublas::zero_matrix<Resource>(2, pb_.nResources());
  }

  // Setup the demands for each Node
  for (Node n : pb_.nodes()) {
    Index i = node_mapping[n.id].id;
    Index p = m_[n].id;
    if (subparts_[p].size() == 1) continue;
    for (Index j = 0; j < pb_.nResources(); ++j) {
      problems[p].demands(i, j) = pb_.demands(n.id, j);
    }
  }

  std::vector<Mapping> next_mappings;
  for (Index i = 0; i < nCurrentParts(); ++i) {
    if (problems[i].nNodes() == 0 || subparts_[i].size() == 1) {
      next_mappings.emplace_back();
      continue;
    }
    // Totals for the subproblem
    Index n_parts_lft = (subparts_[i].size() + 1) / 2;

    for (Index a = 0; a < subparts_[i].size(); ++a) {
      Index p = subparts_[i][a].id;
      Index ind = a >= n_parts_lft;
      for (Index j = 0; j < pb_.nResources(); ++j) {
        problems[i].capacities(ind, j) += pb_.capacities(p, j);
      }
    }

    std::vector<Resource> sub_demands = problems[i].getTotalDemands();
    std::vector<Resource> sub_capacities = problems[i].getTotalCapacities();

    // Setup the capacities, with some margin depending on the size of the subproblem
    // Size 2 ---> no margin
    for (Index j = 0; j < pb_.nResources(); ++j) {
      assert (sub_capacities[j] >= sub_demands[j]);
      // Don't take any chance with floating point
      if (sub_capacities[j] == sub_demands[j] || sub_demands[j] == 0) continue;

      double current_ratio = static_cast<double>(sub_capacities[j]) / sub_demands[j];
      double new_ratio = 1.0 + (current_ratio - 1.0) * 2.0 / subparts_[i].size();
      assert (current_ratio > 1.0 && new_ratio > 1.0);
      double scaling = new_ratio / current_ratio;
      problems[i].capacities(0, j) = problems[i].capacities(0, j) * scaling;
    }

    BipartSolver s(problems[i], options);
    s.run();
    if (!s.success()) {
      std::cerr << "No bipartitioning solution found" << std::endl;
      exit(1);
    }
    next_mappings.emplace_back(s.solution());
  }

  Mapping next_mapping(pb_.nNodes());
  std::vector<std::vector<Part> > next_subparts;

  for (Index i = 0; i < nCurrentParts(); ++i) {
    if (subparts_[i].size() == 1) {
      for (Node n : part_to_nodes[i]) {
        next_mapping[n] = Part(next_subparts.size());
      }
      next_subparts.emplace_back(subparts_[i]);
    }
    else {
      for (Node n : part_to_nodes[i]) {
        next_mapping[n] = Part(next_subparts.size() + next_mappings[i][node_mapping[n.id]].id);
      }
      auto middle = subparts_[i].begin() + (subparts_[i].size() + 1) / 2;
      next_subparts.emplace_back(subparts_[i].begin(), middle);
      next_subparts.emplace_back(middle, subparts_[i].end());
    }
  }

  m_ = next_mapping;
  subparts_ = next_subparts;
}

void BisectionState::checkConsistency() const {
  assert (m_.nNodes() == pb_.nNodes());
}

Mapping solve(const Problem &pb, const SolverOptions &options) {
  BisectionState state(pb);
  state.run(options);

  return state.mapping();
}

} // End namespace minipart


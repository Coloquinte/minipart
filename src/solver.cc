// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#include "solver.h"
#include "bipart_solver.h"

#include <unordered_map>
#include <bitset>
#include <iostream>
#include <fstream>
#include "io.h"

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
  // Bisection tree
  typedef std::vector<std::vector<Index> > LevelToLastLevelParts;
  std::vector<LevelToLastLevelParts> bisection_tree_;
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
  // TODO: initial check that we have a feasible solution with a global bin packing
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
  bisection_tree_.emplace_back();

  // Just setup the bisection (illegal solution)
  for (Index i = 0; i < nCurrentParts(); ++i) {
    Index parts_ind = next_subparts.size();
    if (subparts_[i].size() > 1) {
      auto middle = subparts_[i].begin() + (subparts_[i].size() + 1) / 2;
      next_subparts.emplace_back(subparts_[i].begin(), middle);
      next_subparts.emplace_back(middle, subparts_[i].end());
      next_part_to_nodes.emplace_back(part_to_nodes_[i]);
      next_part_to_nodes.emplace_back();

      bisection_tree_.back().push_back(std::vector<Index>({parts_ind, parts_ind + 1}));
    }
    else {
      next_subparts.emplace_back(subparts_[i]);
      next_part_to_nodes.emplace_back(part_to_nodes_[i]);

      bisection_tree_.back().push_back(std::vector<Index>({parts_ind}));
    }
  }

  subparts_ = next_subparts;
  part_to_nodes_ = next_part_to_nodes;
}

void BisectionState::legalizeBisection() {
  // Run a bin-packing pass on all bisections to get a legal placement
  // If we can't make it at this stage, attempt to question previous bisection levels

  std::vector<char> illegalSubpart(nCurrentParts(), true);
  for (unsigned i = bisection_tree_.size(); i > 0; --i) {
    bool failed_parts = 0;
    const std::vector<std::vector<Index> > &level = bisection_tree_[i-1];
    for (unsigned j = 0; j < level.size(); ++j) {
      // Gather all the parts in it
      std::vector<Index> parts = level[j];
      for (unsigned k = i; k < bisection_tree_.size(); ++k) {
        std::vector<Index> next_parts;
        for (Index p : parts) {
          next_parts.insert(next_parts.end(), bisection_tree_[k][p].begin(), bisection_tree_[k][p].end());
        }
        std::swap(next_parts, parts);
      }

      // Check whether the solution is already legal
      bool illegal = false;
      for (Index p : parts) {
        illegal |= illegalSubpart[p];
      }
      if (!illegal) continue;

      // Legalize them
      bool success = legalizePartitions(parts);
      for (Index p : parts) {
        illegalSubpart[p] = !success;
      }
      if (!success) {
        failed_parts += parts.size();
      }
    }

    if (failed_parts == 0) {
      break;
    }
    else if (options_.verbosity >= 2) {
      std::cout << "Unable to reuse the partitioning at level " << i << ": " << failed_parts << " illegal parts. Trying the next level" << std::endl;
    }
  }

  for (bool illegal : illegalSubpart) {
    if (!illegal) continue;
    std::cerr << "No legal solution found" << std::endl;
    exit(1);
  }
}

void BisectionState::optimizeBisection() {
  for (const std::vector<Index> bisection : bisection_tree_.back()) {
    assert (bisection.size() <= 2 && !bisection.empty());
    if (bisection.size() == 1) continue;
    redoBisection(bisection[0], bisection[1]);
  }
}

void BisectionState::refineBisection() {
  // Up one level
  if (bisection_tree_.size() <= 1) return;

  std::size_t ind = bisection_tree_.size() - 2;
  std::vector<std::vector<Index> > last_level = bisection_tree_.back();
  std::vector<std::vector<Index> > upper_level = bisection_tree_[ind];

  for (const std::vector<Index> bisection : upper_level) {
    // We get 2 x 1 or 2 partitions
    assert (bisection.size() <= 2 && !bisection.empty());
    if (bisection.size() == 1) continue;
    Index u1 = bisection[0];
    Index u2 = bisection[1];
    for (Index p1 : last_level[u1]) {
      for (Index p2 : last_level[u2]) {
        redoBisection(p1, p2);
      }
    }
  }
}

void BisectionState::redoBisection(Index i, Index j) {
  assert (i != j);

  // Extract the problem itself
  BisectionProblem bisection = getBisectionProblem(i, j);

  // Run bipartitioning
  BipartSolver s(bisection.problem, options_);
  // Use initial legal solution
  s.addInitialSolution(bisection.mapping);
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

bool BisectionState::legalizePartitions(const std::vector<Index> &subparts) {
  checkConsistency();
  typedef std::pair<Node, float> P;

  // Get the resources
  Matrix<Resource> capacities = boost::numeric::ublas::zero_matrix<Resource>(subparts.size(), pb_.nResources());
  for (Index i = 0; i < subparts.size(); ++i) {
    for (Part p : subparts_[subparts[i]]) {
      for (unsigned j = 0; j < pb_.nResources(); ++j) {
        capacities(i, j) += pb_.capacities(p.id, j);
      }
    }
  }

  // Normalize each resource
  std::vector<float> factors;
  for (std::size_t j = 0; j < capacities.size2(); ++j) {
    Resource res = 0;
    for (std::size_t i = 0; i < capacities.size1(); ++i) {
      res += capacities(i, j);
    }
    factors.push_back(1.0f / (std::max((float) res, 1.0f)));
  }

  // Get the nodes; sort them by size
  std::vector<P> wNodes;
  for (Index sp : subparts) {
    for (Node n : part_to_nodes_[sp]) {
      float dem_avg = 0.0f;
      for (std::size_t i = 0; i < pb_.nResources(); ++i) {
        dem_avg += pb_.demands(n.id, i) * factors[i];
      }
      wNodes.emplace_back(n, dem_avg);
    }
  }
  std::sort (wNodes.begin(), wNodes.end(),
      [](P a, P b) { return a.second > b.second; });
  std::vector<Node> nodes;
  for (P p : wNodes) nodes.push_back(p.first);

  // Run a greedy bin packing
  std::vector<std::vector<Node> > result(subparts.size());
  Matrix<Resource> left = capacities;
  for (Node n : nodes) {
    std::vector<float> partUsage(capacities.size1(), 0.0);
    for (unsigned j = 0; j < pb_.nResources(); ++j) {
      Resource usage = 0;
      for (unsigned i = 0; i < capacities.size1(); ++i) {
        usage += (capacities(i, j) - left(i, j));
      }
      partUsage[j] += factors[j] * usage;
    }

    // Pick the most "balanced" partition
    Index bestPart = -1;
    float bestUsage = 0.0;
    for (unsigned i = 0; i < capacities.size1(); ++i) {
      bool canPlace = true;
      for (unsigned j = 0; j < pb_.nResources(); ++j) {
        if (pb_.demands(n.id, j) > left(i, j)) canPlace = false;
      }
      if (!canPlace) continue;
      if (bestPart == (Index) -1 || partUsage[i] < bestUsage) {
        bestUsage = partUsage[i];
        bestPart = i;
      }
    }
    if (bestPart == (Index) -1) return false;
    for (unsigned j = 0; j < pb_.nResources(); ++j) {
       left(bestPart, j) -= pb_.demands(n.id, j);
       result[bestPart].push_back(n);
    }
  }

  // Update the solution
  for (Index i = 0; i < subparts.size(); ++i) {
    part_to_nodes_[subparts[i]] = result[i];
  }

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

  std::unordered_map<Index, Index> globalToInternal;
  for (Index i = 0; i < ret.nodes.size(); ++i) {
    globalToInternal[ret.nodes[i].id] = i;
  } 

  // Setup the demands
  ret.problem.demands = Matrix<Resource>(ret.nodes.size(), pb_.nResources());
  for (Index k = 0; k < ret.nodes.size(); ++k) {
    Node n = ret.nodes[k];
    for (Index l = 0; l < pb_.nResources(); ++l) {
      ret.problem.demands(k, l) = pb_.demands(n.id, l);
    }
  }

  // Compute the full capacities
  Matrix<Resource> capacities = boost::numeric::ublas::zero_matrix<Resource>(2, pb_.nResources());
  for (Part p : subparts_[i]) {
    for (Index k = 0; k < pb_.nResources(); ++k) {
      capacities(0, k) += pb_.capacities(p.id, k);
    }
  }
  for (Part p : subparts_[j]) {
    for (Index k = 0; k < pb_.nResources(); ++k) {
      capacities(1, k) += pb_.capacities(p.id, k);
    }
  }

  // Compute the usage for the starting solution
  Matrix<Resource> startingUsage = boost::numeric::ublas::zero_matrix<Resource>(2, pb_.nResources());
  for (Index k = 0; k < ret.nodes.size(); ++k) {
    Node n = ret.nodes[k];
    for (Index l = 0; l < pb_.nResources(); ++l) {
      startingUsage(ret.mapping[Node(k)].id, l) += pb_.demands(n.id, l);
    }
  }
  for (Index k = 0; k < 2; ++k) {
    for (Index l = 0; l < pb_.nResources(); ++l) {
      assert (startingUsage(k, l) <= capacities(k, l));
    }
  }

  // Adjust the capacities: leave some margin for later bisections to succeed
  Matrix<Resource> margins = capacities - startingUsage;
  ret.problem.capacities = Matrix<Resource>(2, pb_.nResources());
  std::size_t num_subparts[2] = { subparts_[i].size(), subparts_[j].size() };
  for (int k = 0; k < 2; ++k) {
    if (num_subparts[k] == 1) {
      // No bisection left: no tweaking required
      for (Index l = 0; l < pb_.nResources(); ++l) {
        ret.problem.capacities(k, l) = capacities(k, l);
      }
    }
    else {
      // Scale the margin
      double scaling = 1.0 / num_subparts[k]; // TODO: think and ponder about it, or benchmark stuff
      for (Index l = 0; l < pb_.nResources(); ++l) {
        ret.problem.capacities(k, l) = startingUsage(k, l) + scaling * margins(k, l);
      }
    }
  }

  // Setup the hypergraph
  HypergraphBuilder b(ret.nodes.size());

  // Gather the hyperedges
  std::unordered_map<Index, std::vector<Node> > edge_to_nodes;
  for (Index k = 0; k < ret.nodes.size(); ++k) {
    Node n = ret.nodes[k];
    Node o(k);
    assert (globalToInternal.at(n.id) == o.id);
    for (Edge e : pb_.hypergraph.edges(n)) {
      edge_to_nodes[e.id].push_back(o);
    }
  }

  // Translate the hyperedges
  for (auto const &p : edge_to_nodes) {
    Edge e(p.first);
    if (p.second.size() == 1) continue;

    // If the cost function is just the cut, only use edges that are not cut
    if (!options_.soed_objective && p.second.size() != pb_.hypergraph.nodes(e).size()) continue;

    b.addEdge(p.second.begin(), p.second.end(), pb_.hypergraph.weight(e));
  }

  // Translate the 2-edges
  for (Index k = 0; k < ret.nodes.size(); ++k) {
    Node n = ret.nodes[k];
    Node o(k);
    for (auto e : pb_.hypergraph.edges2(n)) {
      Node t = e.target;
      if (t.id >= n.id) continue;
      auto it = globalToInternal.find(t.id);
      if (it == globalToInternal.end()) continue;
      b.addEdge2(o, Node(it->second), e.weight);
    }
  }

  b.vectorize();
  ret.problem.hypergraph = b;

  return ret;
}

void BisectionState::checkConsistency() const {
  assert (part_to_nodes_.size() == subparts_.size());
  std::size_t sz = 0;
  for (const auto &nodes : part_to_nodes_) {
    sz += nodes.size();
  }
  assert (sz == pb_.nNodes());
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


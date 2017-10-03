// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#include "inc_bipart.h"

#include <iostream>
#include <iomanip>

namespace minipart {

void place(IncBipart &inc, std::minstd_rand &rgen) {
  std::bernoulli_distribution dist;
  for (auto n : inc.nodes()) {
    if (dist(rgen)) inc.move(n);
  }

  std::vector<Node> nodes (inc.nodes().begin(), inc.nodes().end());
  std::shuffle(nodes.begin(), nodes.end(), rgen);
  for (auto n : nodes) {
    bool m = inc.mapping(n);
    if (inc.overflow(m)) inc.move(n);
  }

  inc.checkConsistency();
}

void optimize(IncBipart &inc, std::minstd_rand &rgen) {
  std::vector<Node > nodes (inc.nodes().begin(), inc.nodes().end());
  std::shuffle(nodes.begin(), nodes.end(), rgen);

  std::vector<Node > active_set, zero_gain_set;
  std::vector<int> num_zero_gain_moves(inc.nNodes(), 0);
  const int max_zero_gain_moves = 2;

  auto tryPush = [&](Node o) {
    if (inc.gain(o) > 0) {
      active_set.push_back(o);
    }
    else if (inc.gain(o) == 0 && ++num_zero_gain_moves[o.id] <= max_zero_gain_moves) {
      zero_gain_set.push_back(o);
    }
  };

  for (auto n : nodes) {
    tryPush(n);
  }

  while (!active_set.empty() || !zero_gain_set.empty()) {
    std::vector<Node > &pick = active_set.empty() ? zero_gain_set : active_set;
    Node n = pick.back();
    pick.pop_back();
    if (inc.gain(n) < 0) continue;
    inc.tryMove(n, [&](Node o, int w) { tryPush(o); });
  }
}

void incrementCutEdges(std::vector<int> &cutEdges, const IncBipart &pb) {
  assert (cutEdges.size() == pb.nEdges());
  for (std::size_t i = 0; i < pb.nEdges(); ++i) {
    if (!pb.cut(Edge(i))) continue;
    ++cutEdges[i];
  }
}

std::size_t getCutUnderCount (const std::vector<int> &cutEdges, int count) {
    std::size_t nEdges = 0;
    for (int nbCut : cutEdges) {
      if (nbCut <= count) {
        ++nEdges;
      }
    }
    return nEdges;
}

void reportCutProportion (const std::vector<int> &cutEdges, int nStarts) {
  std::cout << "\nCut edges: " << std::endl;
  for (double percentage = 0.5; percentage <= 50; percentage *= 2) {
    int maxNbCut = (percentage * 0.01) * nStarts;
    int nb = getCutUnderCount(cutEdges, maxNbCut);
    std::cout << "<= " << percentage << "%: " << 100.0 * nb / cutEdges.size() << "%" << std::endl;
  }
}

std::pair <double, double> computeAvgAndDev(const std::vector<int> &costs) {
  double avg = 0.0;
  double sqavg = 0.0;
  for (auto d : costs) {
    double c = d;
    avg += c;
    sqavg += c * c;
  }
  avg /= costs.size();
  sqavg /= costs.size();
  double std_dev = 100.0 * std::sqrt (sqavg - avg * avg) / avg;
  return std::make_pair(avg, std_dev);
}

void solve(const Problem &pb) {
  std::minstd_rand rgen;
  IncBipart inc(pb);

  const int n_iter = 500;
  std::vector<int> init_cost, final_cost;
  std::vector<int> cut_edges(inc.nEdges(), 0);

  for (int i = 0; i < n_iter; ++i) {
    place(inc, rgen);
    if (!inc.legal()) continue;
    init_cost.push_back(inc.cost());
    optimize(inc, rgen);
    final_cost.push_back(inc.cost());
    incrementCutEdges(cut_edges, inc);
  }

  std::cout << n_iter<< " iterations, ";
  std::cout << std::fixed << std::setw(10) << std::setprecision(2);
  std::cout << 100.0 * (n_iter- init_cost.size()) / n_iter << "% failure, " << std::endl;

  auto init_summary = computeAvgAndDev(init_cost);
  auto final_summary = computeAvgAndDev(final_cost);
  std::cout << "Init:  average " << init_summary.first  << ", deviation " << init_summary.second  << "%" << std::endl;
  std::cout << "Final: average " << final_summary.first  << ", deviation " << final_summary.second  << "%" << std::endl;

  reportCutProportion(cut_edges, n_iter);
}

} // End namespace minipart


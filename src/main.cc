// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#include "hypergraph.h"
#include "io.h"
#include "inc_bipart.h"

#include <boost/program_options.hpp>

#include <fstream>
#include <iostream>
#include <iomanip>

using namespace minipart;
namespace po = boost::program_options;

po::options_description getOptions() {
  po::options_description desc("Options");
  desc.add_options()("help,h", "print this help");

  desc.add_options()("hmetis", po::value<std::string>()->required(),
      "filename for the .hgr file describing the hypergraph");

  desc.add_options()("margin", po::value<double>()->default_value(5.0),
      "margin compared to balanced partitioning (%)")(
      "num_parts", po::value<std::size_t>()->default_value(2u),
      "number of partitions");

  desc.add_options()("dump-hmetis", po::value<std::string>(),
      "dump a .hgr file for debug purposes");
  desc.add_options()("stats", "print graph statistics");

  return desc;
}

void printHelp(const po::options_description &desc) {
  std::cout << "Minipart help:" << std::endl;
  std::cout << desc << std::endl;
}

po::variables_map processOptions(const po::options_description &desc, int argc, char **argv) {
  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
  } catch (po::error &e) {
    std::cerr << "Error parsing command line arguments:" << std::endl;
    std::cerr << e.what() << std::endl;
    printHelp(desc);
    exit(1);
  }
  return vm;
}

po::variables_map parseArguments(int argc, char **argv) {
  po::options_description desc = getOptions();
  po::variables_map vm = processOptions(desc, argc, argv);

  if (vm.count("help")) {
    printHelp(desc);
    exit(0);
  }
  assert(vm.count("margin"));
  assert(vm.count("num_parts"));

  return vm;
}

Problem parseGraph(const po::variables_map &vm) {
  std::ifstream hf(vm["hmetis"].as<std::string>());
  return readHMetis(hf);
}

void exportGraph(const po::variables_map &vm, const Problem &pb) {
  std::ofstream hf(vm["dump-hmetis"].as<std::string>());
  writeHMetis(pb, hf);
}

void setupCapacities(const po::variables_map &vm, Problem &pb) {
  std::vector<int> totals(pb.demands.size2(), 0);
  for (std::size_t i = 0; i < pb.demands.size1(); ++i) {
    for (std::size_t j = 0; j < pb.demands.size2(); ++j) {
      totals[j] += pb.demands(i, j);
    }
  }

  std::size_t num_parts = vm["num_parts"].as<std::size_t>();
  double margin = vm["margin"].as<double>();
  double factor = (margin * 0.01 + 1.0) / num_parts;
  pb.capacities.resize(num_parts, pb.demands.size2());
  for (std::size_t i = 0; i < num_parts; ++i) {
    for (std::size_t j = 0; j < pb.demands.size2(); ++j) {
      pb.capacities(i, j) = totals[j] * factor;
    }
  }
}

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

int main(int argc, char **argv) {
  po::variables_map vm = parseArguments(argc, argv);
  Problem pb = parseGraph(vm);
  setupCapacities(vm, pb);

  if (vm.count("dump-hmetis")) exportGraph(vm, pb);

  if (vm.count("stats")) reportStats(pb.hypergraph, std::cout);

  solve(pb);

  return 0;
}




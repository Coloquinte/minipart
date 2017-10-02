// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#include "hypergraph.h"
#include "io.h"
#include "inc_bipart.h"

#include <boost/program_options.hpp>

#include <fstream>
#include <iostream>

using namespace minipart;
namespace po = boost::program_options;

typedef Problem<unsigned, int, int> PB;

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

PB parseGraph(const po::variables_map &vm) {
  std::ifstream hf(vm["hmetis"].as<std::string>());
  return readHMetis<unsigned, int, int>(hf);
}

void exportGraph(const po::variables_map &vm, const PB &pb) {
  std::ofstream hf(vm["dump-hmetis"].as<std::string>());
  writeHMetis(pb, hf);
}

void setupCapacities(const po::variables_map &vm, PB &pb) {
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

void place(IncBipart<unsigned, int, int> &inc, std::minstd_rand &rgen) {
  std::bernoulli_distribution dist;
  for (auto n : inc.nodes()) {
    if (dist(rgen)) inc.move(n);
  }

  std::vector<Node<unsigned> > nodes (inc.nodes().begin(), inc.nodes().end());
  std::shuffle(nodes.begin(), nodes.end(), rgen);
  for (auto n : nodes) {
    bool m = inc.mapping(n);
    if (inc.overflow(m)) inc.move(n);
  }

  inc.checkConsistency();
}

void optimize(IncBipart<unsigned, int, int> &inc, std::minstd_rand &rgen) {
  std::vector<Node<unsigned> > nodes (inc.nodes().begin(), inc.nodes().end());
  std::shuffle(nodes.begin(), nodes.end(), rgen);

  for (auto n : nodes) {
    if (inc.gain(n) > 0) {
      inc.tryMove(n);
    }
  }

  std::vector<Node<unsigned> > active_set, zero_gain_set;
  std::vector<int> num_zero_gain_moves(inc.nNodes(), 0);
  const int max_zero_gain_moves = 2;

  auto tryPush = [&](Node<unsigned> o) {
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
    Node<unsigned> n;
    if (!active_set.empty()) {
      n = active_set.back();
      active_set.pop_back();
    }
    else {
      n = zero_gain_set.back();
      zero_gain_set.pop_back();
    }
    if (inc.gain(n) < 0) continue;
    inc.tryMove(n, [&](Node<unsigned> o, int w) { tryPush(o); });
  }
}

void solve(const PB &pb) {
  std::minstd_rand rgen;
  const int nTry = 100;
  int nSuccess = 0;

  IncBipart<unsigned, int, int> inc(pb);
  for (int i = 0; i < nTry; ++i) {
    place(inc, rgen);
    if (!inc.legal()) continue;
    ++nSuccess;
    std::cout << "Before " << inc.cost() << std::endl;
    optimize(inc, rgen);
    std::cout << "After " << inc.cost() << std::endl;
  }

  std::cout << nSuccess << " successes out of " << nTry << " placements" << std::endl;
}

int main(int argc, char **argv) {
  po::variables_map vm = parseArguments(argc, argv);
  PB pb = parseGraph(vm);
  setupCapacities(vm, pb);

  if (vm.count("dump-hmetis")) exportGraph(vm, pb);

  if (vm.count("stats")) reportStats(pb.hypergraph, std::cout);

  solve(pb);

  return 0;
}




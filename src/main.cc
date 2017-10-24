// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#include "solver.h"
#include "io.h"

#include <boost/program_options.hpp>

#include <fstream>
#include <iostream>
#include <iomanip>

using namespace minipart;
namespace po = boost::program_options;

void check_margin(double p) {
  if (p < 0.0 || !std::isfinite(p)) {
    throw std::runtime_error("Option --margin must be non-negative");
  }
}

void check_parts(std::size_t p) {
  if (p < 2) {
    throw std::runtime_error("Option --parts must be 2 or greater");
  }
  else if (p > 2) {
    throw std::runtime_error("Only bipartitioning is supported at this point; other values for option --parts are not available");
  }
}

void check_starts(std::size_t p) {
  if (p == 0) {
    throw std::runtime_error("Option --starts must be 1 or greater");
  }
}

void check_cycles(std::size_t p) {
  if (p == 0) {
    throw std::runtime_error("Option --v-cycles must be 1 or greater");
  }
}

void check_threads(std::size_t p) {
  if (p == 0) {
    throw std::runtime_error("Option --threads must be 1 or greater");
  }
}

po::options_description getOptions() {
  po::options_description desc("Minipart options");
  desc.add_options()("help,h", "print this help");

  desc.add_options()("hmetis,g", po::value<std::string>(),
      ".hgr file describing the hypergraph");

  desc.add_options()("output,o", po::value<std::string>(),
      "partitioning result file");

  desc.add_options()("margin", po::value<double>()->default_value(5.0)
      ->notifier(check_margin),
      "margin compared to balanced partitioning (%)");

  desc.add_options()("parts", po::value<std::size_t>()->default_value(2u)
      ->notifier(check_parts),
      "number of partitions");

  desc.add_options()("starts", po::value<std::size_t>()->default_value(32)
      ->notifier(check_starts),
      "number of starting points");

  desc.add_options()("v-cycles", po::value<std::size_t>()->default_value(2)
      ->notifier(check_cycles),
      "number of coarsening-uncoarsening cycles");

  desc.add_options()("threads,j", po::value<std::size_t>()->default_value(2)
      ->notifier(check_threads),
      "number of concurrent threads");

  desc.add_options()("seed", po::value<std::size_t>()->default_value(1),
      "random generator seed");

  desc.add_options()("stats", "print problem statistics");

  return desc;
}

po::options_description getHiddenOptions() {
  po::options_description desc("Hidden options");

  desc.add_options()("dump-hmetis", po::value<std::string>(),
      "dump a .hgr file for debug purposes");

  return desc;
}

po::variables_map parseArguments(int argc, char **argv) {
  po::options_description desc = getOptions();
  po::options_description hidden = getHiddenOptions();

  po::options_description all_options;
  all_options.add(desc).add(hidden);
  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, all_options), vm);
    po::notify(vm);
  } catch (po::error &e) {
    std::cerr << "Error parsing command line arguments: ";
    std::cerr << e.what() << std::endl << std::endl;
    std::cout << desc << std::endl;
    exit(1);
  }
  catch (std::runtime_error &e) {
    std::cerr << "Error parsing command line arguments: ";
    std::cerr << e.what() << std::endl << std::endl;
    exit(1);
  }

  assert(vm.count("margin"));
  assert(vm.count("parts"));
  assert(vm.count("starts"));
  assert(vm.count("seed"));
  assert(vm.count("v-cycles"));

  if (vm.count("help")) {
    std::cout << desc << std::endl;
    exit(0);
  }
  else if (!vm.count("hmetis") || vm["hmetis"].as<std::string>().empty()) {
    std::cout << "Missing input file" << std::endl;
    std::cout << desc << std::endl;
    exit(1);
  }
  else if(vm["parts"].as<std::size_t>() != 2) {
    std::cout << "Only bipartitioning is supported" << std::endl;
    exit(1);
  }

  return vm;
}

Problem parseGraph(const po::variables_map &vm) {
  std::ifstream hf(vm["hmetis"].as<std::string>());
  return readHMetis(hf);
}

void setupCapacities(const po::variables_map &vm, Problem &pb) {
  std::vector<int> totals(pb.demands.size2(), 0);
  for (std::size_t i = 0; i < pb.demands.size1(); ++i) {
    for (std::size_t j = 0; j < pb.demands.size2(); ++j) {
      totals[j] += pb.demands(i, j);
    }
  }

  std::size_t num_parts = vm["parts"].as<std::size_t>();
  double margin = vm["margin"].as<double>();
  double factor = (margin * 0.01 + 1.0) / num_parts;
  pb.capacities.resize(num_parts, pb.demands.size2());
  for (std::size_t i = 0; i < num_parts; ++i) {
    for (std::size_t j = 0; j < pb.demands.size2(); ++j) {
      pb.capacities(i, j) = totals[j] * factor;
    }
  }
}

void reportInputs(const po::variables_map &vm, Problem &pb) {
  if (vm.count("dump-hmetis")) {
    std::ofstream hf(vm["dump-hmetis"].as<std::string>());
    writeHMetis(pb, hf);
  }
  if (vm.count("stats")) reportStats(pb, std::cout);
}

void writeResult(const po::variables_map &vm, const std::vector<Mapping> &mappings) {
  if (!mappings.empty() && vm.count("output")) {
    std::ofstream of(vm["output"].as<std::string>());
    for (std::size_t i = 0; i < mappings.front().nNodes(); ++i) {
      of << (std::size_t) mappings.front()[Node(i)].id << std::endl;
    }
  }
}

void check_solutions(const Problem &pb, std::vector<Mapping> &pool) {
  // TODO: add checkers
}

void sort_solutions(const Problem &pb, std::vector<Mapping> &pool) {
  typedef std::pair<std::int64_t, Mapping> CM;
  std::vector<CM> cost_to_mapping;
  for (const Mapping & m : pool) {
    // TODO: get rid of bipart-specific stuff
    cost_to_mapping.emplace_back(computeBipartCost(pb.hypergraph, m), m);
  }
  std::sort (cost_to_mapping.begin(), cost_to_mapping.end(), [](const CM &a, const CM &b) { return a.first < b.first; });
  pool.clear();
  for (CM &c : cost_to_mapping) {
    pool.emplace_back(c.second);
  }
}

int main(int argc, char **argv) {
  po::variables_map vm = parseArguments(argc, argv);
  Problem pb = parseGraph(vm);
  setupCapacities(vm, pb);

  SolverOptions opt;
  opt.n_starts  = vm["starts"].as<std::size_t>();
  opt.n_cycles  = vm["v-cycles"].as<std::size_t>();
  opt.n_threads = vm["threads"].as<std::size_t>();
  opt.seed      = vm["seed"].as<std::size_t>();

  reportInputs(vm, pb);

  std::vector<Mapping> mappings = solve(pb, opt);

  check_solutions(pb, mappings);
  sort_solutions(pb, mappings);

  reportResults(pb, mappings, std::cout);

  writeResult(vm, mappings);

  return 0;
}




// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#include "solver.h"
#include "io.h"

#include <boost/program_options.hpp>

#include <fstream>
#include <iostream>
#include <iomanip>

using namespace minipart;
namespace po = boost::program_options;

po::options_description getOptions() {
  po::options_description desc("Minipart options");
  desc.add_options()("help,h", "print this help");

  desc.add_options()("hmetis", po::value<std::string>(),
      "filename for the .hgr file describing the hypergraph");

  desc.add_options()("margin", po::value<double>()->default_value(5.0),
      "margin compared to balanced partitioning (%)");

  desc.add_options()("parts", po::value<std::size_t>()->default_value(2u),
      "number of partitions");

  desc.add_options()("starts", po::value<int>()->default_value(32),
      "number of starting points");

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

  if (vm.count("help")) {
    std::cout << desc << std::endl;
    exit(0);
  }
  else if (!vm.count("hmetis") || vm["hmetis"].as<std::string>().empty()) {
    std::cout << "Missing input file" << std::endl;
    std::cout << desc << std::endl;
    exit(1);
  }

  assert(vm.count("margin"));
  assert(vm.count("parts"));
  assert(vm.count("starts"));
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
  if (vm.count("dump-hmetis")) exportGraph(vm, pb);
  if (vm.count("stats")) reportStats(pb, std::cout);
}

int main(int argc, char **argv) {
  po::variables_map vm = parseArguments(argc, argv);
  Problem pb = parseGraph(vm);
  setupCapacities(vm, pb);
  reportInputs(vm, pb);

  std::vector<Mapping> mappings = solve(pb, vm["starts"].as<int>());
  reportResults(pb, mappings, std::cout);

  return 0;
}




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
}

void check_objective(std::string obj) {
  if (obj != "cut" && obj!= "degree") {
    throw std::runtime_error("Option --objective must be \"cut\" ot \"degree\"");
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

  desc.add_options()("objective", po::value<std::string>()->default_value("cut")
      ->notifier(check_objective),
      "objective cost function");

  desc.add_options()("starts", po::value<std::size_t>()->default_value(64)
      ->notifier(check_starts),
      "number of starting points");

  desc.add_options()("v-cycles", po::value<std::size_t>()->default_value(5)
      ->notifier(check_cycles),
      "number of coarsening-uncoarsening cycles");

  desc.add_options()("threads,j", po::value<std::size_t>()->default_value(4)
      ->notifier(check_threads),
      "number of concurrent threads");

  desc.add_options()("seed", po::value<std::size_t>()->default_value(0),
      "random generator seed");

  desc.add_options()("verbose,v", po::value<std::size_t>()->default_value(1),
      "verbosity level");

  return desc;
}

po::options_description getHiddenOptions() {
  po::options_description desc("Hidden options");

  desc.add_options()("stats", "print problem statistics");
  desc.add_options()("export-hmetis", po::value<std::string>(), "export to hmetis format");

  return desc;
}

po::variables_map parse_arguments(int argc, char **argv) {
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
  assert(vm.count("verbose"));
  assert(vm.count("v-cycles"));

  if (vm.count("help")) {
    std::cout << desc << std::endl;
    exit(0);
  }
  if (!vm.count("hmetis") || vm["hmetis"].as<std::string>().empty()) {
    std::cout << "Missing input file" << std::endl;
    std::cout << desc << std::endl;
    exit(1);
  }
  if (vm["parts"].as<std::size_t>() > 64) {
    std::cout << "The maximum number of partitions is 64" << std::endl;
    exit(1);
  }

  return vm;
}

Problem parse_graph(const po::variables_map &vm) {
  std::ifstream hf(vm["hmetis"].as<std::string>());
  Problem ret = readHMetis(hf);
  if (vm.count("export-hmetis")) {
    std::ofstream of(vm["export-hmetis"].as<std::string>());
    writeHMetis(ret, of);
  }
  return ret;
}

void setup_capacities(const po::variables_map &vm, Problem &pb) {
  std::vector<int> totals = pb.getTotalDemands();

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

void report_inputs(const po::variables_map &vm, Problem &pb) {
  if (vm.count("stats")) reportStats(pb, std::cout);
  pb.check_consistency();
}

void write_solution(const po::variables_map &vm, const Mapping &mapping) {
  if (vm.count("output")) {
    std::ofstream of(vm["output"].as<std::string>());
    for (std::size_t i = 0; i < mapping.nNodes(); ++i) {
      of << (std::size_t) mapping[Node(i)].id << std::endl;
    }
  }
}

void preprocess(Problem &pb) {
  bool error = false;
  for (Node n : pb.nodes()) {
    bool foundValidPlacement = false;
    for (std::size_t i = 0; i < pb.nParts(); ++i) {
      bool valid = true;
      for (std::size_t j = 0; j < pb.nResources(); ++j) {
        if (pb.demands(n.id, j) > pb.capacities(i, j)) valid = false;
      }
      if (valid) foundValidPlacement = true;
    }
    if (!foundValidPlacement) {
      std::cerr << "Node " << n.id << " cannot be placed on any partition" << std::endl;
      error = true;
    }
  }
  if (error) exit(1);
}

void check_solution(const Problem &pb, const Mapping &sol) {
  if (!pb.is_legal(sol)) {
    std::cerr << "Illegal solution returned" << std::endl;
    exit(1);
  }
}

SolverOptions get_options(const po::variables_map &vm) {
  SolverOptions opt;
  opt.n_starts  = vm["starts"].as<std::size_t>();
  opt.n_cycles  = vm["v-cycles"].as<std::size_t>();
  opt.n_threads = vm["threads"].as<std::size_t>();
  opt.seed      = vm["seed"].as<std::size_t>();
  opt.soed_objective = vm["objective"].as<std::string>() == "degree";
  opt.verbosity = vm["verbose"].as<std::size_t>();

  return opt;
}

int main(int argc, char **argv) {
  po::variables_map vm = parse_arguments(argc, argv);
  Problem pb = parse_graph(vm);
  setup_capacities(vm, pb);
  report_inputs(vm, pb);
  preprocess(pb);

  SolverOptions opt = get_options(vm);
  Mapping mapping = solve(pb, opt);
  write_solution(vm, mapping);
  check_solution(pb, mapping);

  if (opt.verbosity >= 1) {
    std::cout << computeCostCut(pb.hypergraph, mapping);
    if (pb.nParts() > 2) std::cout << " " << computeCostDegree(pb.hypergraph, mapping);
    std::cout << std::endl;
  }
  return 0;
}




// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#include "io.h"

#include <ostream>
#include <istream>
#include <iomanip>
#include <map>

namespace minipart {

namespace {
std::stringstream getHMetisLine(std::istream &s) {
  // Get rid of comment lines and empty lines
  std::string tmp;
  do {
    if (!s.good()) throw std::runtime_error("Not enough lines");
    std::getline(s, tmp);
  } while(tmp.empty() || tmp[0] == '%');

  return std::stringstream(tmp);
}

Hypergraph readHMetisGraph(std::istream &s, std::size_t nNodes, std::size_t nEdges, bool hasEdgeWeights) {
  std::stringstream ss;
  HypergraphBuilder h(nNodes);

  std::vector<Node> nodes;
  for (std::size_t i = 0; i < nEdges; ++i) {
    ss = getHMetisLine(s);

    Weight w = 1;
    if (hasEdgeWeights) ss >> w;

    while (ss) {
        Index n;
        ss >> n;
        if (ss.fail()) continue;
        if (n > nNodes) throw std::runtime_error("Parsed pin index is outside the specified number of nodes");
        if (n == 0) throw std::runtime_error("Parsed pin index cannot be 0");
        nodes.emplace_back(n-1);
    }
    if (nodes.empty()) throw std::runtime_error("No node on the line");

    h.addEdge(nodes.begin(), nodes.end(), w);

    nodes.clear();
  }

  return h;
}

Matrix<Resource> readHMetisResources(std::istream &s, std::size_t nNodes, bool hasNodeWeights) {
  std::stringstream ss;
  Matrix<Resource> r(nNodes, 1);
  if (hasNodeWeights) {
    std::vector<Resource> resources;
    for (std::size_t i = 0; i < nNodes; ++i) {
      ss = getHMetisLine(s);
      while (ss) {
        Resource w;
        ss >> w;
        if (ss.fail()) continue;
        resources.push_back(w);
      }
      if (i == 0) r = Matrix<Resource>(nNodes, resources.size());
      if (resources.size() != r.size2()) throw std::runtime_error("Inconsistent number of node weights");
      for (std::size_t j = 0; j < r.size2(); ++j) {
        r(i, j) = resources[j];
      }
      resources.clear();
    }
  }
  else {
    for (std::size_t i = 0; i < nNodes; ++i) {
      r(i, 0) = 1;
    }
  }

  return r;
}
} // End anonymous namespace

Problem readHMetis(std::istream &s) {
  Problem ret;

  std::size_t nNodes, nEdges, params;

  std::stringstream ss = getHMetisLine(s);
  ss >> nEdges >> nNodes;
  if (ss.fail()) throw std::runtime_error("Invalid first line");
  ss >> params;
  if (ss.fail()) params = 0;

  if (params != 0 && params != 1 && params != 10 && params != 11) throw std::runtime_error("Invalid parameter value");

  bool hasEdgeWeights = (params == 11) || (params == 1);
  bool hasNodeWeights = (params == 11) || (params == 10);

  ret.hypergraph = readHMetisGraph(s, nNodes, nEdges, hasEdgeWeights);
  ret.demands = readHMetisResources(s, nNodes, hasNodeWeights);
  return ret;
}

void writeHMetis(const Problem &pb, std::ostream &s) {
  s << "% HGR (hMetis) file generated by Minipart" << std::endl;

  const Hypergraph &h = pb.hypergraph;
  const Matrix<Resource> &r = pb.demands;

  s << "% " << h.nNodes() << " nodes, " << h.nEdges() << " edges" << std::endl;
  if (r.size2() > 1) {
    s << "% Multiple (" << r.size2() << ") node weights: this file might not be readable by other tools" << std::endl;
  }
  s << "%" << std::endl;

  s << h.nEdges() << " " << h.nNodes();
  s << " 11" << std::endl; // Edge and node weights

  for (auto e : h.edges()) {
    s << h.weight(e);
    for (auto n : h.nodes(e)) {
      s << " " << n.id+1;
    }
    s << std::endl;
  }

  for (std::size_t i = 0; i < r.size1(); ++i) {
    for (std::size_t j = 0; j < r.size2(); ++j) {
      s << r(i, j) << " ";
    }
    s << std::endl;
  }
}

void reportHelper (
    std::ostream &s
  , const std::map<std::size_t, std::size_t> &degreeToCount
  , const std::map<std::size_t, std::size_t> &degreeToWeight
  , bool reportWeight = false
  ) {
  std::size_t totCnt = 0;
  std::size_t totPins = 0;
  std::size_t totWeight = 0;

  for (auto p : degreeToCount) {
    totCnt += p.second;
    totPins += p.first * p.second;
  }
  for (auto p : degreeToWeight) {
    totWeight += p.second;
  }

  for (std::size_t i = 0u; i < 8u; ++i) {
    if (degreeToCount.count(i) == 0) continue;
    s << i << ",\t";
    s << 100.0 * degreeToCount.at(i) / totCnt << "%,\t";
    s << 100.0 * i * degreeToCount.at(i) / totPins << "%";
    if (reportWeight) s << ",\t" << 100.0 * degreeToWeight.at(i) / totWeight << "%";
    s << std::endl;
  }

  std::size_t maxIndex = degreeToCount.rbegin()->first;
  for (std::size_t i = 8u; i < maxIndex; i *= 2) {
    std::size_t cnt = 0;
    std::size_t pins = 0;
    std::size_t wgt = 0;
    for (std::size_t j = i; j < 2 * i; ++j) {
      if (degreeToCount.count(j) == 0) continue;
      cnt += degreeToCount .at(j);
      pins += degreeToCount.at(j) * j;
      wgt += degreeToWeight.at(j);
    }
    s << i << "-" << 2 * i - 1 << ",\t";
    s << 100.0 * cnt / totCnt << "%,\t";
    s << 100.0 * pins / totPins << "%";
    if (reportWeight) s << ",\t" << 100.0 * wgt / totWeight << "%";
    s << std::endl;
  }

  s << std::endl;
}

void reportStats(const Problem &pb, std::ostream &s) {
  const Hypergraph &h = pb.hypergraph;
  std::map<std::size_t, std::size_t> degreeToCount;
  std::map<std::size_t, std::size_t> degreeToWeight;

  std::size_t sumDegrees = 0;
  bool reportWeight = false;

  for (auto e : h.edges()) {
    std::size_t degree = h.nodes(e).size();
    std::size_t weight = h.weight(e);
    ++degreeToCount[degree];
    degreeToWeight[degree] += weight;
    sumDegrees += degree;
    if (weight != 1) reportWeight = true;
  }
  if (degreeToCount.empty()) return;

  s << std::endl;
  s << h.nNodes() << " nodes, " << h.nEdges() << " edges" << std::endl;
  s << std::fixed << std::setw(3) << std::setprecision(1);
  s << static_cast<double>(sumDegrees) / h.nEdges() << " avg edge degree" << std::endl;
  s << static_cast<double>(sumDegrees) / h.nNodes() << " avg node degree" << std::endl;
  s << std::endl;

  s << "Degree,\tEdges,\tPins";
  if (reportWeight) s << ",\tWeight";
  s << std::endl;

  reportHelper(s, degreeToCount, degreeToWeight, reportWeight);

  degreeToCount.clear();
  for (auto n : h.nodes()) {
    std::size_t degree = h.edges(n).size();
    ++degreeToCount[degree];
  }
  degreeToWeight = degreeToCount;

  s << "Degree,\tNodes,\tPins";
  s << std::endl;

  reportHelper(s, degreeToCount, degreeToWeight);
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

bool isEdgeCut(const Hypergraph &h, const Mapping &m, Edge e) {
  // Used as a bitset
  std::uint8_t used = 0;
  for (auto n : h.nodes(e)) {
    uint8_t pos = m[n].id;
    assert (pos <= 1);
    used |= (pos + 1);
  }
  assert (used <= 3);
  // Cut if both bits are true
  return used == 3;
}

std::int64_t computeBipartCost(const Hypergraph &h, const Mapping &m) {
  std::int64_t ret = 0;
  for (auto e : h.edges()) {
    if (isEdgeCut(h, m, e)) {
      ret += h.weight(e);
    }
  }
  return ret;
}

std::vector<int> getCutCounts(const Hypergraph &h, const std::vector<Mapping> &mappings) {
  std::vector<int> cut_counts(h.nEdges(), 0);
  for (Edge e : h.edges()) {
    for (const Mapping &m : mappings) {
      cut_counts[e.id] += isEdgeCut(h, m, e);
    }
  }
  return cut_counts;
}

std::size_t getCutUnderCount (const std::vector<int> &cut_counts, int max_count) {
    std::size_t num = 0;
    for (int cut_count : cut_counts) {
      num += (cut_count <= max_count);
    }
    return num;
}

void reportResults(const Problem &pb, const std::vector<Mapping> &mappings, std::ostream &s) {
  std::cout << mappings.size() << " feasible placements, ";

  std::vector<int> costs;
  for (const Mapping &m : mappings) {
    costs.push_back(computeBipartCost(pb.hypergraph, m));
  }

  auto summary = computeAvgAndDev(costs);
  s << std::fixed << std::setw(10) << std::setprecision(2);
  s << "Cost: average " << summary.first  << ", deviation " << summary.second  << "%" << std::endl;

  std::vector<int> cut_counts = getCutCounts(pb.hypergraph, mappings);
  s << "\nCut edges: " << std::endl;
  for (double percentage = 0.5; percentage <= 50; percentage *= 2) {
    int max_count = (percentage * 0.01) * mappings.size();
    int nb = getCutUnderCount(cut_counts, max_count);
    s << "<= " << percentage << "%: " << 100.0 * nb / cut_counts.size() << "%" << std::endl;
  }

}

}  // End namespace minipart


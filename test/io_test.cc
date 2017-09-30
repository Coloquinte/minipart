
#define BOOST_TEST_MODULE PARSE_TEST
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "hypergraph.h"
#include "io.h"
#include "rand_gen.h"

#include <sstream>
#include <iostream>

using namespace minipart;

typedef Hypergraph<unsigned, int> H;

BOOST_AUTO_TEST_SUITE(io)

BOOST_AUTO_TEST_CASE(hypergraph) {
  unsigned nNodes = 10;
  unsigned nEdges = 10;
  unsigned nResources = 2;
  std::minstd_rand rgen;

  Problem<unsigned, int, int> pb;
  H::Builder b(nNodes);
  addRandomEdges(b, nEdges, 3.0, 10.0, rgen, 1u);
  pb.hypergraph = b;
  pb.hypergraph.checkConsistency();
  pb.demands = getRandomArray<int>(nNodes, nResources, 10, rgen);

  std::stringstream s;
  writeHMetis(pb, s);
  Problem<unsigned, int, int> copied = readHMetis<unsigned, int, int>(s);
  BOOST_CHECK_EQUAL (pb.hypergraph.nNodes(), copied.hypergraph.nNodes());
  BOOST_CHECK_EQUAL (pb.hypergraph.nEdges(), copied.hypergraph.nEdges());
  BOOST_CHECK_EQUAL (pb.hypergraph.nPins(),  copied.hypergraph.nPins());
  for (auto e : pb.hypergraph.edges()) {
    BOOST_CHECK_EQUAL (pb.hypergraph.nodes(e),  copied.hypergraph.nodes(e));
    BOOST_CHECK_EQUAL (pb.hypergraph.weight(e),  copied.hypergraph.weight(e));
  }

  BOOST_CHECK_EQUAL (pb.demands.size1(), copied.demands.size1());
  BOOST_CHECK_EQUAL (pb.demands.size2(), copied.demands.size2());
  for (unsigned i = 0; i < nNodes; ++i) {
    for (unsigned j = 0; j < nResources; ++j) {
      BOOST_CHECK_EQUAL (pb.demands(i, j), copied.demands(i, j));
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()


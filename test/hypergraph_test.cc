
#define BOOST_TEST_MODULE HYPERGRAPH_TEST
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "hypergraph.h"

using namespace minipart;

BOOST_AUTO_TEST_SUITE(hypergraph)

BOOST_AUTO_TEST_CASE(degree2) {
  typedef Hypergraph<unsigned, int> H;
  H::Builder b;
  
  std::vector<H::N> nodes;
  for (unsigned i = 0; i < 5; ++i) {
    nodes.push_back(b.addNode());
  }

  std::vector<H::E> edges;
  for (unsigned i = 0; i < 5; ++i) {
    for (unsigned j = i+1; j < 5; ++j) {
      H::E e = b.addEdge({nodes[i], nodes[j]}, edges.size());
      edges.push_back(e);
    }
  }

  H g = b;
  g.checkConsistency();

  BOOST_CHECK_EQUAL (g.nEdges(), 10u);
  BOOST_CHECK_EQUAL (g.nNodes(), 5u);
  BOOST_CHECK (edges.size() == g.nEdges());
  BOOST_CHECK (nodes.size() == g.nNodes());
  for (auto e : edges) {
    BOOST_CHECK_EQUAL (g.weight(e), e.id);
    BOOST_CHECK_EQUAL (g.nodes(e).size(), 2lu);
  }
  for (auto n : g.nodes()) {
    BOOST_CHECK_EQUAL (g.edges(n).size(), 4lu);
  }
}

BOOST_AUTO_TEST_CASE(degree3) {
  typedef Hypergraph<unsigned, int> H;
  H::Builder b;
  
  std::vector<H::N> nodes;
  for (unsigned i = 0; i < 5; ++i) {
    nodes.push_back(b.addNode());
  }

  std::vector<H::E> edges;
  for (unsigned i = 0; i < 5; ++i) {
    for (unsigned j = i+1; j < 5; ++j) {
      for (unsigned k = j+1; k < 5; ++k) {
        H::E e = b.addEdge({nodes[i], nodes[j], nodes[k]}, edges.size());
        edges.push_back(e);
      }
    }
  }

  H g = b;
  g.checkConsistency();

  BOOST_CHECK_EQUAL (g.nEdges(), 10u);
  BOOST_CHECK_EQUAL (g.nNodes(), 5u);
  for (auto e : edges) {
    BOOST_CHECK_EQUAL (g.weight(e), e.id);
    BOOST_CHECK_EQUAL (g.nodes(e).size(), 3lu);
  }
  for (auto n : g.nodes()) {
    BOOST_CHECK_EQUAL (g.edges(n).size(), 6lu);
  }
}

BOOST_AUTO_TEST_SUITE_END()


// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#pragma once

#include <cstddef>
#include <cstdint>
#include <boost/range/iterator_range_core.hpp>
#include <boost/numeric/ublas/matrix.hpp>

namespace minipart {

template <typename T>
using Matrix = boost::numeric::ublas::matrix<T>;

// Typesafe handles: Node/Edge/Part

template <typename Index>
struct Node {
  Index id;

  Node () : id(-1) {}
  explicit Node (Index i) : id(i) {}
  bool operator==(Node o) const { return id == o.id; }
  bool operator!=(Node o) const { return id != o.id; }
  bool operator<(Node o) const { return id < o.id; }

  template <typename S> operator Node<S>() { return Node<S>(id); }
};

template <typename Index>
struct Part {
  Index id;

  Part () : id(-1) {}
  explicit Part (Index i) : id(i) {}
  bool operator==(Part o) const { return id == o.id; }
  bool operator!=(Part o) const { return id != o.id; }
  bool operator<(Part o) const { return id < o.id; }

  template <typename S> operator Part<S>() { return Part<S>(id); }
};

template <typename Index>
struct Edge {
  Index id;

  Edge () : id(-1) {}
  explicit Edge (Index i) : id(i) {}
  bool operator==(Edge o) const { return id == o.id; }
  bool operator!=(Edge o) const { return id != o.id; }
  bool operator<(Edge o) const { return id < o.id; }
};

// Node/Edge iterators
template<class T>
struct ObjIt : T, boost::iterator_facade<ObjIt<T>, T, boost::forward_traversal_tag> {
  ObjIt(T i) : T(i) {}
  ObjIt<T>& operator++() { ++this->id; return *this; }
  T& operator*() { return *this; }
};

template<class T> using Range = boost::iterator_range<ObjIt<T> >;

// All Pin iterators (Node/Reader/Driver, Edge/Input/Output)
template <class T> using SliceIt = typename std::vector<T>::const_iterator;  
template <class T> using Slice = boost::iterator_range<SliceIt<T> >;  

template <class Index>
class Mapping {
 public:
  typedef Part<Index> P;

  Mapping(std::size_t s=0) : m_(s, P(0)) {}
  std::size_t nNodes() const { return m_.size(); }

  template <typename I>
  P  operator[](Node<I> n) const { return m_[n.id]; }
  template <typename I>
  P& operator[](Node<I> n)       { return m_[n.id]; }

 private:
  std::vector<P> m_;
};

}  // End namespace minipart

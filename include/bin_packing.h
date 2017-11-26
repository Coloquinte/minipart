
#pragma once

#include "common.h"

namespace minipart {
class BinPackingSolver {
 public:
  BinPackingSolver(const Matrix<Resource> &demands, const Matrix<Resource> &capacities)
    : demands_(demands)
    , capacities_(capacities) {
    capacity_left_ = capacities_;
  }

  bool success() const { return m_.nNodes() == demands_.size1(); }
  const Mapping &mapping() const { return m_; }

  void run(std::minstd_rand &rgen) {
    std::vector<Node> nodes;
    for (Index i = 0; i < demands_.size1(); ++i) {
      nodes.emplace_back(i);
    }
    std::shuffle(nodes.begin(), nodes.end(), rgen);
    std::stable_sort(nodes.begin(), nodes.end(),
      [&](Node a, Node b) {
        // TODO: handle several resources
        return demands_(a.id, 0) > demands_(b.id, 0);
      }
    );

    m_ = Mapping(demands_.size1());
    for (Node n : nodes) {
      bool ok = false;
      for (Index i = 0; !ok && i < capacities_.size1(); ++i) {
        ok = tryPlace(n, i);
      }
      if (!ok) {
        m_ = Mapping();
        break;
      }
    }
  }

 private:
  bool tryPlace(Node n, Index i) {
    bool canPlace = true;
    for (std::size_t j = 0; j < demands_.size2(); ++j) {
      if (demands_(n.id, j) > capacity_left_(i, j)) {
        canPlace = false;
      }
    }
    if (!canPlace) return false;
    for (std::size_t j = 0; j < demands_.size2(); ++j) {
      capacity_left_(i, j) -= demands_(n.id, j);
    }
    m_[n] = Part(i);
    return true;
  }

 private:
  const Matrix<Resource> &demands_;
  const Matrix<Resource> &capacities_;
  Matrix<Resource> capacity_left_;
  Mapping m_;
};

}


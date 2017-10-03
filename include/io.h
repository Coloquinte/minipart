// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#pragma once

#include "hypergraph.h"

#include <iosfwd>

namespace minipart {

Problem readHMetis(std::istream &s);
void writeHMetis(const Problem &pb, std::ostream &s);

void reportStats(const Hypergraph &h, std::ostream &s);

}  // End namespace minipart


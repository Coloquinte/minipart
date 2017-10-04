// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#pragma once

#include "solver.h"

#include <iosfwd>

namespace minipart {

Problem readHMetis(std::istream &s);
void writeHMetis(const Problem &pb, std::ostream &s);

void reportStats(const Problem &pb, std::ostream &s);
void reportResults(const Problem &pb, const std::vector<Mapping> &m, std::ostream &s);

}  // End namespace minipart


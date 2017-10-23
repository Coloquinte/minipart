// Copyright (C) 2017 Gabriel Gouvine - All Rights Reserved

#pragma once

#include "inc_bipart.h"

namespace minipart {

// Placement passes
void random_placement_pass(IncBipart &inc, std::minstd_rand &rgen);

// Local search passes
void greedy_pass(IncBipart &inc, std::minstd_rand &rgen, int passes=3);
void positive_gain_pass(IncBipart &inc, std::minstd_rand &rgen);
}


#include "multiplier_mcmc_move.h"
#include "beagle_tree_likelihood.h"
#include <cmath>
#include <iostream>
#include <utility>

using namespace bpp;

namespace sts { namespace online {

struct Proposal
{
    double value;
    double hastings_ratio;
};

Proposal pos_real_multiplier(const double, const double, const double, const double, smc::rng*);

}}

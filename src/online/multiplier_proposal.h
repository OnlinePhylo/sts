#include "multiplier_mcmc_move.h"
//#include "beagle_tree_likelihood.h"
#include <cmath>
#include <iostream>
#include <utility>

namespace sts { namespace online {

struct Proposal
{
    double value;
    double forwardDensity;  // For SMC
    double hastingsRatio;   // For MCMC
};

Proposal positive_real_multiplier(const double, const double, const double, const double, smc::rng*);

}}

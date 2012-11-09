#ifndef STS_MOVES_EXPONENTIAL_BRANCH_LENGTH_PROPOSER_HPP
#define STS_MOVES_EXPONENTIAL_BRANCH_LENGTH_PROPOSER_HPP

#include <cassert>
#include <cmath>
#include <utility>
#include "gsl/gsl_randist.h"
#include "smctc.hh"

#include "sts/particle/state.hpp"
#include "sts/moves/base_branch_length_proposer.hpp"

namespace sts
{
namespace moves
{

/// \class exponential_branch_length_proposer
/// \brief Propose branch lengths from an exponential distribution.
class exponential_branch_length_proposer : public base_branch_length_proposer
{
public:
    /// Instantiate a new BL proposer where branch lengths are drawn from an exponential distribution with mean
    /// \c mean.
    /// \param mean Mean of exponential distribution
    explicit exponential_branch_length_proposer(double mean) : mean(mean) {};
    double log_proposal_density(double);

    /// Mean of exponential distribution
    double mean;
protected:
    double propose_bl(smc::rng *rng);
};

/// Propose a branch length

/// \returns A pair consisting of: (branch_length, likelihood)
double exponential_branch_length_proposer::propose_bl(smc::rng *rng)
{
    return rng->Exponential(this->mean);
}

double exponential_branch_length_proposer::log_proposal_density(double d)
{
    return -(std::log(mean) + (d / mean));
}

} // namespace moves
} // namespace sts

#endif // STS_MOVES_EXPONENTIAL_BRANCH_LENGTH_PROPOSER_HPP

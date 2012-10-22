#ifndef STS_MOVES_EXPONENTIAL_BRANCH_LENGTH_PROPOSAL_HPP
#define STS_MOVES_EXPONENTIAL_BRANCH_LENGTH_PROPOSAL_HPP

#include <cassert>
#include <cmath>
#include <utility>
#include "gsl/gsl_randist.h"
#include "smctc.hh"

#include "sts/particle/phylo_particle.hpp"
#include "sts/moves/branch_length_proposer.hpp"

namespace sts
{
namespace moves
{

/// \class exponential_branch_length_proposer
/// \brief Propose branch lengths from an exponential distribution.
class exponential_branch_length_proposer : public branch_length_proposer
{
public:
    /// Instantiate a new BL proposer where branch lengths are drawn from an exponential distribution with mean
    /// \c mean.
    /// \param mean Mean of exponential distribution
    explicit exponential_branch_length_proposer(double mean) : mean(mean) {};
    double log_proposal_density(double);
    branch_length_proposer::branch_lengths propose(particle::particle, smc::rng *);
protected:
    /// Mean of exponential distribution
    double mean;
    double propose_bl(smc::rng *rng);
};

/// Propose branch lengths on \c part.
branch_length_proposer::branch_lengths exponential_branch_length_proposer::propose(particle::particle part, smc::rng *rng)
{
    // TODO: different BLs for child1 and child2
    double d1 = propose_bl(rng); // , d2 = propose_bl(rng)
    return branch_length_proposer::branch_lengths(d1, d1);
}

/// Propose a branch length

/// \returns A pair consisting of: (branch_length, likelihood)
double exponential_branch_length_proposer::propose_bl(smc::rng *rng)
{
    return rng->Exponential(this->mean);
}

double exponential_branch_length_proposer::log_proposal_density(double d)
{
    return std::log(gsl_ran_exponential_pdf(d, this->mean));
}

} // namespace moves
} // namespace sts

#endif // STS_MOVES_EXPONENTIAL_BRANCH_LENGTH_PROPOSAL_HPP

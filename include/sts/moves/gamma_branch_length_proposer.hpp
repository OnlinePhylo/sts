#ifndef STS_MOVES_GAMMA_BRANCH_LENGTH_PROPOSER_HPP
#define STS_MOVES_GAMMA_BRANCH_LENGTH_PROPOSER_HPP

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

/// \class gamma_branch_length_proposer
/// \brief Propose branch lengths from an gamma distribution.
class gamma_branch_length_proposer : public branch_length_proposer
{
public:
    /// Instantiate a new BL proposer where branch lengths are drawn from an gamma distribution with mean
    /// \c mean.
    /// \param mean Mean of gamma distribution
    explicit gamma_branch_length_proposer(double mean) : mean(mean) {};
    double log_proposal_density(double);
    branch_length_proposer::branch_lengths propose(particle::particle, smc::rng *);

    /// Mean of gamma distribution
    double mean;
protected:
    double propose_bl(smc::rng *rng);
};

/// Propose branch lengths on \c part.
branch_length_proposer::branch_lengths gamma_branch_length_proposer::propose(particle::particle part, smc::rng *rng)
{
    // TODO: different BLs for child1 and child2
    double d1 = propose_bl(rng); // , d2 = propose_bl(rng)
    return branch_length_proposer::branch_lengths(d1, d1);
}

/// Propose a branch length

/// \returns A pair consisting of: (branch_length, likelihood)
double gamma_branch_length_proposer::propose_bl(smc::rng *rng)
{
    return rng->Gamma(2.0, this->mean);
}

double gamma_branch_length_proposer::log_proposal_density(double d)
{
    double gamma_p = gsl_ran_gamma_pdf(d, 2.0, this->mean);
    return std::log(gamma_p);
}

} // namespace moves
} // namespace sts

#endif // STS_MOVES_GAMMA_BRANCH_LENGTH_PROPOSER_HPP

#ifndef STS_MOVES_DELTA_BRANCH_LENGTH_PROPOSER_HPP
#define STS_MOVES_DELTA_BRANCH_LENGTH_PROPOSER_HPP

#include <cassert>
#include <cmath>
#include <utility>
#include <limits>
#include "smctc.hh"

#include "sts/particle/phylo_particle.hpp"
#include "sts/moves/branch_length_proposer.hpp"

namespace sts
{
namespace moves
{

/// \class delta_branch_length_proposer
/// \brief "Propose" branch lengths from a delta distribution.
class delta_branch_length_proposer : public branch_length_proposer
{
public:
    /// Instantiate a new BL proposer where branch lengths are drawn from an delta distribution with mean
    /// \c mean.
    /// \param mean Mean of delta distribution
    explicit delta_branch_length_proposer(double mean) : mean(mean) {};
    double log_proposal_density(double);
    branch_length_proposer::branch_lengths propose(particle::particle, smc::rng *);

    /// Mean of delta distribution
    double mean;
protected:
    double propose_bl(smc::rng *rng);
};

/// Propose branch lengths on \c part.
branch_length_proposer::branch_lengths delta_branch_length_proposer::propose(particle::particle part, smc::rng *rng)
{
    // TODO: different BLs for child1 and child2
    double d1 = propose_bl(rng); // , d2 = propose_bl(rng)
    return branch_length_proposer::branch_lengths(d1, d1);
}

/// Propose a branch length

/// \returns A pair consisting of: (branch_length, likelihood)
double delta_branch_length_proposer::propose_bl(smc::rng *rng)
{
    return this->mean;
}

double delta_branch_length_proposer::log_proposal_density(double d)
{
    if(d == this->mean) return 0;
    else return -std::numeric_limits<double>::infinity();
}

} // namespace moves
} // namespace sts

#endif // STS_MOVES_DELTA_BRANCH_LENGTH_PROPOSER_HPP

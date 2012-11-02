#ifndef STS_MOVES_DELTA_BRANCH_LENGTH_PROPOSER_HPP
#define STS_MOVES_DELTA_BRANCH_LENGTH_PROPOSER_HPP

#include <cassert>
#include <cmath>
#include <utility>
#include <limits>
#include "smctc.hh"

#include "sts/particle/phylo_particle.hpp"
#include "sts/moves/base_branch_length_proposer.hpp"

namespace sts
{
namespace moves
{

/// \class delta_branch_length_proposer
/// \brief "Propose" branch lengths from a delta distribution.
class delta_branch_length_proposer : public base_branch_length_proposer
{
public:
    /// Instantiate a new BL proposer where branch lengths are drawn from an delta distribution with mean
    /// \c mean.
    /// \param mean Mean of delta distribution
    explicit delta_branch_length_proposer(double mean) : mean(mean) {};
    double log_proposal_density(double) override;

    /// Mean of delta distribution
    double mean;
protected:
    double propose_bl(smc::rng *rng) override;
};

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

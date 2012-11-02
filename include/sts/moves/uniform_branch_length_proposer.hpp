#ifndef STS_MOVES_UNIFORM_BRANCH_LENGTH_PROPOSER_HPP
#define STS_MOVES_UNIFORM_BRANCH_LENGTH_PROPOSER_HPP

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

/// \class uniform_branch_length_proposer
/// \brief Propose branch lengths from a uniform distribution.
class uniform_branch_length_proposer : public base_branch_length_proposer
{
public:
    /// Instantiate a new BL proposer where branch lengths are drawn from an uniform distribution on [0, 2 \c mean].
    /// \param mean Mean of uniform distribution
    explicit uniform_branch_length_proposer(double mean) : mean(mean) {};
    double log_proposal_density(double) override;

    /// Mean of uniform distribution
    double mean;
protected:
    double propose_bl(smc::rng *rng) override;
};
/// Propose a branch length

/// \returns A pair consisting of: (branch_length, likelihood)
double uniform_branch_length_proposer::propose_bl(smc::rng *rng)
{
    return rng->Uniform(0, 2 * this->mean);
}

double uniform_branch_length_proposer::log_proposal_density(double d)
{
    if(0 <= d && d >= 2 * this->mean) return -log(2*this->mean); // log( 1/(2 mu) )
    else return -std::numeric_limits<double>::infinity();
}

} // namespace moves
} // namespace sts

#endif // STS_MOVES_UNIFORM_BRANCH_LENGTH_PROPOSER_HPP

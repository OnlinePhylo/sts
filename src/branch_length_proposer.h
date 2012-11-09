#ifndef STS_MOVES_BRANCH_LENGTH_PROPOSER_HPP
#define STS_MOVES_BRANCH_LENGTH_PROPOSER_HPP

#include "smctc.hh"

#include "state.h"

namespace sts
{
namespace moves
{

/// \class branch_length_proposal
/// \brief Abstract class
class Branch_length_proposer
{
public:
    virtual ~Branch_length_proposer() {};

    /// Propose branch lengths on \c node.

    /// \param part Phylo node to operate on. Child edges of \c node must be initialized.
    /// <b>This function changes child edge branch lengths.</b>
    /// \param rng Random number generator
    /// \returns The log-likelihood of the proposal
    virtual double operator()(particle::Particle part, smc::rng* rng) = 0;

    /// Prior density for proposal with branch-length d.
    /// \param d Branch length
    /// \returns Log probability of branch with length \c d
    virtual double log_proposal_density(double d) = 0;
};

} // namespace moves
} // namespace sts

#endif // STS_MOVES_BRANCH_LENGTH_PROPOSER_HPP

#ifndef STS_MOVES_BRANCH_LENGTH_PROPOSER_H
#define STS_MOVES_BRANCH_LENGTH_PROPOSER_H

#include "particle.h"
#include "smctc.hh"

namespace sts
{
namespace moves
{

/// \class Branch_length_proposer
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
    double operator()(particle::Particle part, smc::rng* rng) { return this->propose_branches(part, rng); };

    /// \param part Phylo node to operate on. Child edges of \c node must be initialized.
    /// <b>This function changes child edge branch lengths.</b>
    /// \param rng Random number generator
    /// \returns The log-likelihood of the proposal
    virtual double propose_branches(particle::Particle part, smc::rng* rng) = 0;

    /// Prior density for proposal with branch-length d.
    /// \param d Branch length
    /// \returns Log probability of branch with length \c d
    virtual double log_proposal_density(double d) = 0;
};

} // namespace moves
} // namespace sts

#endif // STS_MOVES_BRANCH_LENGTH_PROPOSER_H

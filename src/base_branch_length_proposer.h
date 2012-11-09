#ifndef STS_MOVES_BASE_BRANCH_LENGTH_PROPOSER_HPP
#define STS_MOVES_BASE_BRANCH_LENGTH_PROPOSER_HPP

#include <cmath>
#include <utility>
#include "smctc.hh"

#include "branch_length_proposer.h"

namespace sts
{
namespace moves
{
/// \class Base_branch_length_proposer
/// \brief Abstract class
/// Derived classes should implement \c propose_bl and \c log_proposal_density
class Base_branch_length_proposer
{
public:
    /// Convenience type - pair of branches
    typedef std::pair<double, double> Branch_lengths;

    /// Propose branch lengths on \c node.

    /// \param part Phylo node to operate on. Child edges of \c node must be initialized.
    /// <b>This function changes child edge branch lengths.</b>
    /// \param rng Random number generator
    /// \returns The log-likelihood of the proposal
    double operator()(particle::Particle part, smc::rng* rng);

    /// Prior density for proposal with branch-length d.
    /// \param d Branch length
    /// \returns Log probability of branch with length \c d
    virtual double log_proposal_density(double d) = 0;

    /// Propose a pair of branch lengths
    virtual Branch_lengths propose(particle::Particle, smc::rng *);

    virtual ~Base_branch_length_proposer() {};
protected:
    /// Override in subclass
    virtual double propose_bl(smc::rng *rng) = 0;
};

} // namespace moves
} // namespace sts

#endif // STS_MOVES_BASE_BRANCH_LENGTH_PROPOSER_HPP

#ifndef STS_MOVES_BRANCH_LENGTH_PROPOSER_HPP
#define STS_MOVES_BRANCH_LENGTH_PROPOSER_HPP

#include <cassert>
#include <cmath>
#include <utility>
#include "gsl/gsl_randist.h"
#include "smctc.hh"

#include "sts/particle/phylo_particle.hpp"

namespace sts
{
namespace moves
{

/// \class branch_length_proposal
/// \brief Abstract class
class branch_length_proposer
{
public:
    /// Convenience type.
    typedef std::pair<double, double> branch_lengths;

    /// Propose branch lengths on \c node.

    /// \param part Phylo node to operate on. Child edges of \c node must be initialized. <b>This function changes child
    /// edge branch lengths.</b>
    /// \param rng Random number generator
    /// \returns The log-likelihood of the proposal
    virtual double operator()(particle::particle part, smc::rng* rng);

    /// Prior density for proposal with branch-length d.
    /// \param d Branch length
    /// \returns Log probability of branch with length \c d
    virtual double log_proposal_density(double d) = 0;

    /// Propose a pair of branch lengths
    virtual branch_lengths propose(particle::particle, smc::rng *) = 0;
};

// Implementation
double branch_length_proposer::operator()(particle::particle part, smc::rng *rng)
{
    branch_lengths p = propose(part, rng); // This is where the subclassing action happens.
    std::shared_ptr<particle::phylo_node> node = part->node;

    // Children should be initialized
    assert(node->child1);
    assert(node->child2);
    node->child1->length = p.first;
    node->child2->length = p.second;
    return log_proposal_density(p.first) + log_proposal_density(p.second);
}

} // namespace moves
} // namespace sts

#endif // STS_MOVES_BRANCH_LENGTH_PROPOSER_HPP

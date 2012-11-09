#ifndef STS_MOVES_BASE_BRANCH_LENGTH_PROPOSER_HPP
#define STS_MOVES_BASE_BRANCH_LENGTH_PROPOSER_HPP

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

/// \class base_branch_length_proposer
/// \brief Abstract class
/// Derived classes should implement \c propose_bl and \c log_proposal_density
class base_branch_length_proposer
{
public:
    /// Convenience type - pair of branches
    typedef std::pair<double, double> branch_lengths;

    /// Propose branch lengths on \c node.

    /// \param part Phylo node to operate on. Child edges of \c node must be initialized.
    /// <b>This function changes child edge branch lengths.</b>
    /// \param rng Random number generator
    /// \returns The log-likelihood of the proposal
    double operator()(particle::particle part, smc::rng* rng);

    /// Prior density for proposal with branch-length d.
    /// \param d Branch length
    /// \returns Log probability of branch with length \c d
    virtual double log_proposal_density(double d) = 0;

    /// Propose a pair of branch lengths
    virtual branch_lengths propose(particle::particle, smc::rng *);

    virtual ~base_branch_length_proposer() {};
protected:
    /// Override in subclass
    virtual double propose_bl(smc::rng *rng) = 0;
};

// Implementation
double base_branch_length_proposer::operator()(particle::particle part, smc::rng *rng)
{
    branch_lengths p = propose(part, rng); // This is where the subclassing action happens.
    std::shared_ptr<particle::Node> node = part->node;

    // Children should be initialized
    assert(node->child1);
    assert(node->child2);
    node->child1->length = p.first;
    node->child2->length = p.second;
    return log_proposal_density(p.first) + log_proposal_density(p.second);
}

base_branch_length_proposer::branch_lengths base_branch_length_proposer::propose(particle::particle part, smc::rng *rng)
{
    double d1 = propose_bl(rng), d2 = propose_bl(rng);
    return base_branch_length_proposer::branch_lengths(d1, d2);
}

} // namespace moves
} // namespace sts

#endif // STS_MOVES_BASE_BRANCH_LENGTH_PROPOSER_HPP

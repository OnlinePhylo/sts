#ifndef STS_MOVES_BASE_BRANCH_LENGTH_PROPOSER_H
#define STS_MOVES_BASE_BRANCH_LENGTH_PROPOSER_H

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
/// Derived classes should implement \c propose_single_branch_length and \c log_proposal_density
class Base_branch_length_proposer : public Branch_length_proposer
{
public:
    /// Convenience type - pair of branches
    typedef std::pair<double, double> Branch_lengths;

    double propose_branches(particle::Particle part, smc::rng* rng);

    /// Prior density for proposal with branch-length d.
    /// \param d Branch length
    /// \returns Log probability of branch with length \c d
    virtual double log_proposal_density(double d) = 0;

    /// Propose a pair of branch lengths
    virtual Branch_lengths propose(particle::Particle, smc::rng *);

    virtual ~Base_branch_length_proposer() {};
protected:
    /// Override in subclass
    virtual double propose_single_branch_length(smc::rng *rng) = 0;
};

} // namespace moves
} // namespace sts

#endif // STS_MOVES_BASE_BRANCH_LENGTH_PROPOSER_H

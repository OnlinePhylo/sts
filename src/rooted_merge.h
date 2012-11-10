#ifndef STS_MOVES_ROOTED_MERGE_H
#define STS_MOVES_ROOTED_MERGE_H
#include <functional>

#include "forest_likelihood.h"
#include "exponential_branch_length_proposer.h"
#include "smc_move.h"
#include "state.h"

namespace sts
{
namespace moves
{

/// \class Rooted_merge
/// \brief Merge of two nodes, with exponential branch length proposal
class Rooted_merge: public Smc_move
{
public:
    /// Branch length proposal function.
    /// Accepts two parameters: a Node with initialized edges and a random source; returns the log-likelihood.
    typedef std::function<double(particle::Particle, smc::rng*)> Bl_proposal_fn;

    /// Constructor

    /// Initializes with exponential_branch_length_proposal with mean 1.0.
    explicit Rooted_merge(sts::likelihood::Forest_likelihood& log_likelihood) : Smc_move(log_likelihood),
        bl_proposal(Exponential_branch_length_proposer(1.0)) {};

    /// Constructor

    /// \param bl_proposal Source of branch length proposals.
    Rooted_merge(sts::likelihood::Forest_likelihood& log_likelihood,
                 Bl_proposal_fn bl_proposal) : Smc_move(log_likelihood), bl_proposal(bl_proposal) {};

    int do_move(long, smc::particle<particle::Particle>&, smc::rng*) const;

protected:
    /// Branch length proposal generator
    Bl_proposal_fn bl_proposal;
};
} // namespace moves
} // namespace sts

#endif // STS_MOVES_ROOTED_MERGE_H

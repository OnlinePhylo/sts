#ifndef STS_MOVES_ROOTED_MERGE_H
#define STS_MOVES_ROOTED_MERGE_H
#include <functional>

#include "branch_length_proposer.h"
#include "forest_likelihood.h"
#include "exponential_branch_length_proposer.h"
#include "uniform_pair_proposer.h"
#include "smc_move.h"
#include "state.h"

namespace sts
{
namespace moves
{

/// \brief Merge of two nodes, with supplied proposal function.
class Rooted_merge: public Smc_move
{
public:

    /// Node pair proposal function.
    typedef std::function<void(particle::Particle, smc::rng*, particle::Node_ptr&, particle::Node_ptr&, double&, double&)> Pair_proposal_fn;

    /// Constructor

    /// \param log_likelihood  Forest likelihood
    /// \param proposer Source of branch length proposals.
    Rooted_merge(sts::likelihood::Forest_likelihood* log_likelihood,
                 Branch_length_proposer* proposer, Pair_proposal_fn pair_proposal) : Smc_move(log_likelihood), bl_proposer(proposer), pair_proposal(pair_proposal) {};

    int do_move(long, smc::particle<particle::Particle>&, smc::rng*) const;

protected:
    /// Branch length proposal generator
    Branch_length_proposer* bl_proposer;
    /// Node pair proposal generator
    Pair_proposal_fn pair_proposal;
};
} // namespace moves
} // namespace sts

#endif // STS_MOVES_ROOTED_MERGE_H

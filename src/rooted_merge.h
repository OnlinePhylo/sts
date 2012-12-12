#ifndef STS_MOVES_ROOTED_MERGE_H
#define STS_MOVES_ROOTED_MERGE_H
#include <functional>

#include "bl_proposal_fn.h"
#include "forest_likelihood.h"
#include "exponential_branch_length_proposer.h"
#include "uniform_pair_proposer.h"
#include "smc_move.h"
#include "state.h"

namespace sts
{
namespace moves
{

/// \class Rooted_merge
/// \brief Merge of two nodes, with supplied proposal function.
class Rooted_merge: public Smc_move
{
public:

    /// Node pair proposal function.
    typedef std::function<void(particle::Particle, smc::rng*, particle::Node_ptr&, particle::Node_ptr&, double&, double&)> Pair_proposal_fn;

    /// Constructor

    /// Initializes with exponential_branch_length_proposal with mean 1.0.
    explicit Rooted_merge(sts::likelihood::Forest_likelihood& log_likelihood) : Smc_move(log_likelihood),
        bl_proposal(Exponential_branch_length_proposer(1.0)), pair_proposal(Uniform_pair_proposer(log_likelihood)) {};

    /// Constructor

    /// \param bl_proposal Source of branch length proposals.
    Rooted_merge(sts::likelihood::Forest_likelihood& log_likelihood,
                 Bl_proposal_fn bl_proposal, Pair_proposal_fn pair_proposal) : Smc_move(log_likelihood), bl_proposal(bl_proposal), pair_proposal(pair_proposal) {};

    int do_move(long, smc::particle<particle::Particle>&, smc::rng*) const;

protected:
    /// Branch length proposal generator
    Bl_proposal_fn bl_proposal;
    /// Node pair proposal generator
    Pair_proposal_fn pair_proposal;
};
} // namespace moves
} // namespace sts

#endif // STS_MOVES_ROOTED_MERGE_H

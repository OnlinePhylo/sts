#ifndef STS_MOVES_ROOTED_MERGE_H
#define STS_MOVES_ROOTED_MERGE_H
#include <functional>

#include "branch_length_proposer.h"
#include "forest_likelihood.h"
#include "exponential_branch_length_proposer.h"
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
    /// Constructor

    /// \param bl_proposal Source of branch length proposals.
    Rooted_merge(sts::likelihood::Forest_likelihood* log_likelihood,
                 Branch_length_proposer* proposer) : Smc_move(log_likelihood), bl_proposer(proposer) {};

    int do_move(long, smc::particle<particle::Particle>&, smc::rng*) const;

protected:
    /// Branch length proposal generator
    Branch_length_proposer* bl_proposer;
};
} // namespace moves
} // namespace sts

#endif // STS_MOVES_ROOTED_MERGE_H

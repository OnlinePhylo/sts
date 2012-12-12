/// \file child_swap_mcmc_move.h
/// \brief Child_swap_mcmc_move class

#ifndef STS_MOVES_CHILD_SWAP_MCMC_MOVE_H
#define STS_MOVES_CHILD_SWAP_MCMC_MOVE_H

#include "branch_length_proposer.h"
#include "bl_proposal_fn.h"

#include "metropolis_hastings_move.h"

#include "smctc.hh"
#include "state.h"
#include "node.h"

#include <memory>

namespace sts
{
namespace moves
{

/// An MCMC move that replaces one of the two children of this node with an uncoalesced node.
class Child_swap_mcmc_move : public Metropolis_hastings_move
{
public:
    Child_swap_mcmc_move(sts::likelihood::Forest_likelihood& log_likelihood,
                         Bl_proposal_fn* bl_proposer) :
        Metropolis_hastings_move(log_likelihood),
        branch_length_proposer(bl_proposer) {};

    void propose_move(long time, particle::Particle& part, smc::rng* rng) const;
protected:
    Bl_proposal_fn *branch_length_proposer;
};

} // namespace moves
} // namespace sts

#endif // STS_MOVES_CHILD_SWAP_MCMC_MOVE_H

/// \file child_swap_mcmc_move.h
/// \brief Child_swap_mcmc_move class

#ifndef STS_MOVES_CHILD_SWAP_MCMC_MOVE_H
#define STS_MOVES_CHILD_SWAP_MCMC_MOVE_H

#include "mcmc_move.h"

#include "smctc.hh"
#include "state.h"
#include "node.h"

namespace sts
{
namespace moves
{

/// An MCMC move that replaces one of the two children of this node with an uncoalesced node.
class Child_swap_mcmc_move : public Mcmc_move
{
public:
    Child_swap_mcmc_move(sts::likelihood::Forest_likelihood& log_likelihood) : Mcmc_move(log_likelihood) {};

    void propose_move(long time, particle::Particle& part, smc::rng* rng) const;
};

} // namespace moves
} // namespace sts

#endif // STS_MOVES_CHILD_SWAP_MCMC_MOVE_H

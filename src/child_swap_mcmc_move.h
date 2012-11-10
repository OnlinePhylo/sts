/// \file child_swap_mcmc_move.h
/// \brief Child_swap_mcmc_move class

#ifndef STS_MOVES_CHILD_SWAP_MCMC_MOVE_H
#define STS_MOVES_CHILD_SWAP_MCMC_MOVE_H

#include "mcmc_move.h"

#include "smctc.hh"
#include "mcmc_move.h"
#include "state.h"
#include "node.h"

namespace sts
{
namespace moves
{

/// An MCMC move which perturbs branch lengths uniformly from -amount to amount
class Child_swap_mcmc_move : public Mcmc_move
{
public:
    Child_swap_mcmc_move(sts::likelihood::Forest_likelihood& log_likelihood) : Mcmc_move(log_likelihood) {};

    int do_move(long, smc::particle<particle::Particle>&, smc::rng*) const;
};
} // namespace moves
} // namespace sts

#endif // STS_MOVES_CHILD_SWAP_MCMC_MOVE_H
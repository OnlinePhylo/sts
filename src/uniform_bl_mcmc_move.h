/// \file uniform_bl_mcmc_move.h
/// \brief Uniform_bl_mcmc_move class

#ifndef STS_MOVES_UNIFORM_BL_MCMC_MOVE_H
#define STS_MOVES_UNIFORM_BL_MCMC_MOVE_H

#include "smctc.hh"
#include "metropolis_hastings_move.h"
#include "state.h"
#include "node.h"

namespace sts
{
namespace moves
{

/// An MCMC move which perturbs branch lengths uniformly from -amount to amount.
class Uniform_bl_mcmc_move : public Metropolis_hastings_move
{
public:

    Uniform_bl_mcmc_move(sts::likelihood::Forest_likelihood* log_likelihood, double amount) :
        Metropolis_hastings_move(log_likelihood), amount(amount) {};
    Uniform_bl_mcmc_move(sts::likelihood::Forest_likelihood* log_likelihood) :
        Uniform_bl_mcmc_move(log_likelihood, 0.1) {};
    void propose_move(long time, particle::Particle& part, smc::rng* rng) const;

    /// Amount to perturb branch lengths.
    double amount;
};
} // namespace moves
} // namespace sts

#endif // STS_MOVES_UNIFORM_BL_MCMC_MOVE_H

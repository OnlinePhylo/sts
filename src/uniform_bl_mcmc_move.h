/// \file uniform_bl_mcmc_move.h
/// \brief Uniform_bl_mcmc_move class

#ifndef STS_MOVES_UNIFORM_BL_MCMC_MOVE_H
#define STS_MOVES_UNIFORM_BL_MCMC_MOVE_H

#include "smctc.hh"
#include "mcmc_move.h"
#include "state.h"
#include "node.h"

namespace sts
{
namespace moves
{

/// An MCMC move which perturbs branch lengths uniformly from -amount to amount.
class Uniform_bl_mcmc_move : public Mcmc_move
{
public:
    Uniform_bl_mcmc_move(sts::likelihood::Forest_likelihood& log_likelihood) : Mcmc_move(log_likelihood), amount(0.1) {};
    Uniform_bl_mcmc_move(sts::likelihood::Forest_likelihood& log_likelihood, double amount) : Mcmc_move(log_likelihood), amount(amount) {};

    void propose_move(long time, particle::Particle& part, smc::rng* rng) const;

    /// Amount to perturb branch lengths.
    double amount;
};
} // namespace moves
} // namespace sts

#endif // STS_MOVES_UNIFORM_BL_MCMC_MOVE_H

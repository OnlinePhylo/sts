/// \file uniform_bl_mcmc_move.hpp
/// \brief uniform_bl_mcmc_move class

#ifndef STS_MOVES_UNIFORM_BL_MCMC_MOVE_HPP
#define STS_MOVES_UNIFORM_BL_MCMC_MOVE_HPP

#include "smctc.hh"
#include "sts/moves/mcmc_move.hpp"
#include "sts/particle/phylo_particle.hpp"
#include "sts/particle/phylo_node.hpp"

namespace sts
{
namespace moves
{

/// An MCMC move which perturbs branch lengths uniformly from -amount to amount
class uniform_bl_mcmc_move : public mcmc_move
{
public:
    uniform_bl_mcmc_move(sts::likelihood::forest_likelihood& log_likelihood) : mcmc_move(log_likelihood), amount(0.1) {};
    uniform_bl_mcmc_move(sts::likelihood::forest_likelihood& log_likelihood, double amount) : mcmc_move(log_likelihood), amount(amount) {};

    int do_move(long, smc::particle<particle::particle>&, smc::rng*) const;

    /// Amount to perturb branch lengths
    double amount;
};

/// Uniform change to branch lengths for the current node

/// Change the branch lengths for the current node by drawing from a uniform distribution between -amount and amount.
///  \param time  generation number
///  \param from  Source particle
///  \param rng   Random number source
int uniform_bl_mcmc_move::do_move(long time, smc::particle<particle::particle>& from, smc::rng* rng) const
{
    auto calc = log_likelihood.get_calculator();
    particle::particle part = *from.GetValuePointer();
    particle::node cur_node = part->node;
    particle::node new_node = std::make_shared<particle::phylo_node>(calc);
    new_node->child1 = cur_node->child1;
    new_node->child2 = cur_node->child2;

    double cur_ll = log_likelihood(part);
    // Choose an amount to shift the node height uniformly at random.
    double shift = rng->Uniform(-amount, amount);
    // If the shift amount would create a negative node height we will reflect it back to a positive number.
    // This means the proposal density for the reflection zone is double the rest of the area, but the back-proposal
    // probability is also double in the same area so these terms cancel in the Metropolis-Hastings ratio.

    // Now calculate the new node heights - shift both heights for now: ultrametric
    new_node->child1->length = std::abs(new_node->child1->length + shift);
    new_node->child2->length = std::abs(new_node->child2->length + shift);
    part->node = new_node;

    double alpha = exp(log_likelihood(part) - cur_ll);
    if(alpha < 1 && rng->UniformS() > alpha) {
        // Move rejected, restore the original node.
        part->node = cur_node;
        return false;
    }
    // Accept the new state.
    return true;
}

} // namespace moves
} // namespace sts

#endif // STS_MOVES_UNIFORM_BL_MCMC_MOVE_HPP

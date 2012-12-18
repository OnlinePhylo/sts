#include "uniform_bl_mcmc_move.h"
#include "edge.h"

#include <cmath>
#include <cstdlib>
#include <memory>

namespace sts
{
namespace moves
{

/// Uniform change to branch lengths for the current node.
/// Change the branch lengths for the current node by drawing from a uniform distribution between -amount and amount.
///  \param time  generation number
///  \param part  Source particle
///  \param rng   Random number source
double Uniform_bl_mcmc_move::propose_move(long int time, particle::Particle& part, smc::rng* rng) const
{
    // Choose an amount to shift the node height uniformly at random.
    double shift = rng->Uniform(-amount, amount);
    // If the shift amount would create a negative node height we will reflect it back to a positive number.
    // This means the proposal density for the reflection zone is double the rest of the area, but the back-proposal
    // probability is also double in the same area so these terms cancel in the Metropolis-Hastings ratio.

    // Now calculate the new node heights - shift both heights for now: ultrametric
    part->node->child1->length = std::abs(part->node->child1->length + shift);
    part->node->child2->length = std::abs(part->node->child2->length + shift);
    return shift;
}


std::string Uniform_bl_mcmc_move::get_name() const
{
    return "Uniform_bl_mcmc_move";
}

} // namespace moves
} // namespace sts

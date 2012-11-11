#include "mcmc_move.h"

namespace sts
{
namespace moves
{
/// Function call for use with smctc - calls user-defined do_move, tracks result.

/// \param time Generation number
/// \param from Source particle
/// \param rng Random number source
int Mcmc_move::operator()(long time, smc::particle<particle::Particle>& from, smc::rng *rng)
{
    attempted++;

    auto calc = log_likelihood.get_calculator();
    particle::Particle part = *from.GetValuePointer();
    double cur_ll = log_likelihood(part);
    particle::Node_ptr cur_node = part->node;
    particle::Node_ptr new_node = std::make_shared<particle::Node>(*cur_node);

    part->node = new_node;

    // Propose the MCMC move
    propose_move(time, part, rng);
    const double new_ll = log_likelihood(part);

    // Metropolis-Hastings step
    double alpha = exp(new_ll - cur_ll);
    if(alpha < 1 && rng->UniformS() > alpha) {
        // Move rejected, restore the original node.
        part->node = cur_node;
        new_node.reset();
        return false;
    }
    // Accept the new state.
    accepted++;
    return true;
}

} // namespace sts::moves
} // namespace sts

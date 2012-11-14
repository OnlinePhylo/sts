#include "mcmc_move.h"

#include "edge.h"

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
    particle::Node_ptr cur_node = part->node;
    const double cur_ll = log_likelihood(part);
    const double cur_prior = cur_node->child1->prior_log_likelihood + cur_node->child2->prior_log_likelihood;
    particle::Node_ptr new_node = std::make_shared<particle::Node>(*cur_node);

    part->node = new_node;

    // Propose the MCMC move
    propose_move(time, part, rng);
    const double new_ll = log_likelihood(part);
    const double new_prior = new_node->child1->prior_log_likelihood + new_node->child2->prior_log_likelihood;

    // Metropolis-Hastings step
    double alpha = exp(new_ll + new_prior - cur_ll - cur_prior);
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

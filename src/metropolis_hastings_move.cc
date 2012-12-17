#include "metropolis_hastings_move.h"

#include "edge.h"

using namespace sts::particle;

namespace sts
{
namespace moves
{
/// Function call for use with smctc - calls user-defined propose_move, determines whether the move should be accepted.

/// \param time Generation number
/// \param from Source particle
/// \param rng Random number source
/// \returns \c true if the move was accepted.
int Metropolis_hastings_move::operator()(long time, smc::particle<particle::Particle>& from, smc::rng *rng)
{
    attempted++;

    auto calc = log_likelihood->get_calculator();
    particle::Particle cur_part = from.GetValue();
    particle::Particle new_part = std::make_shared<State>();
    const double cur_ll = log_likelihood->calculate_log_likelihood(cur_part);

    // Under a uniform topological prior, the prior for the current node is proportional to the branch length prior.
    const double cur_prior = cur_part->node->edge_prior_log_likelihood();
    particle::Node_ptr new_node = std::make_shared<particle::Node>(*cur_part->node);

    new_part->node = new_node;

    // Propose the MH move
    propose_move(time, new_part, rng);
    const double new_ll = log_likelihood->calculate_log_likelihood(new_part);
    const double new_prior = new_node->edge_prior_log_likelihood();

    // Metropolis-Hastings step
    double alpha = exp(new_ll + new_prior - cur_ll - cur_prior);
    if(alpha < 1 && rng->UniformS() > alpha) {
        // Move rejected - no changes to `from`
        return false;
    }
    // Accept the new state.
    from.SetValue(new_part);
    accepted++;
    return true;
}

/// Probability of this operator being accepted.
double Metropolis_hastings_move::acceptance_probability() const
{
    return (double)attempted / (double)(attempted + accepted);
}

} // namespace sts::moves
} // namespace sts

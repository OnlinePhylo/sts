#include "metropolis_hastings_move.h"

#include "edge.h"

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
    particle::Particle part = *from.GetValuePointer();
    particle::Node_ptr cur_node = part->node;
    const double cur_ll = log_likelihood->calculate_log_likelihood(part);
    
    particle::Particle new_part = std::make_shared<sts::particle::State>();
    new_part->predecessor = part->predecessor;

    // Under a uniform topological prior, the prior for the current node is proportional to the branch length prior.
    const double cur_prior = cur_node->edge_prior_log_likelihood();
    new_part->node = std::make_shared<particle::Node>(*cur_node);
    from.SetValue(new_part);

    // Propose the MH move
    propose_move(time, new_part, rng);
    const double new_ll = log_likelihood->calculate_log_likelihood(new_part);
    const double new_prior = new_part->node->edge_prior_log_likelihood();

    // Metropolis-Hastings step
    double alpha = exp(new_ll + new_prior - cur_ll - cur_prior);
    if(alpha < 1 && rng->UniformS() > alpha) {
        // Move rejected, restore the original node.
	from.SetValue(part);
        new_part->node.reset();
        return false;
    }
    // Accept the new state.
    accepted++;
    return true;
}

/// Probability of this operator being accepted.
double Metropolis_hastings_move::acceptance_probability() const
{
    return (double)attempted / (double)accepted;
}

} // namespace sts::moves
} // namespace sts

#include "rooted_merge.h"
#include "edge.h"
#include "node.h"
#include "util.h"

#include <cassert>
#include <cmath>
#include <functional>
#include <memory>
#include <vector>

namespace sts
{
namespace moves
{
/// Perform a particle move.
/// Merge trees in \c p_from in place, and updates their weights.
/// \f[
/// w_r(s_r) = \frac{\gamma_r(s_r)}{\gamma_{r-1}(s_{r-1})} \frac{\nu^-(s_r \rightarrow s_{r-1})}{\nu^+(s_{r-1} \rightarrow s_r)}.
/// \f]
///\param time The sampler iteration.
///\param p_from The particle to move.
///\param rng  A random number generator.
int Rooted_merge::do_move(long time, smc::particle<particle::Particle>& p_from, smc::rng* rng) const
{
    auto calc = log_likelihood->get_calculator();
    particle::Particle prev_part = p_from.GetValue();
    double prev_ll = log_likelihood->calculate_log_likelihood(prev_part);
    particle::Particle pp = std::make_shared<particle::State>();
    pp->predecessor = p_from.GetValue();
    p_from.SetValue(pp);
    pp->node = std::make_shared<particle::Node>(calc);

    // Select a pair of nodes to merge
    particle::Node_ptr n1, n2;
    double fwd_density, back_density;
    pair_proposal(pp, rng, n1, n2, fwd_density, back_density);

    assert(n1 != nullptr);
    assert(n2 != nullptr);

    // Draw branch lengths.
    pp->node->child1 = std::make_shared<particle::Edge>(n1);
    pp->node->child2 = std::make_shared<particle::Edge>(n2);

    // Because the proposal distribution is uniform on merges, the only part of the proposal distribution we have to
    // worry about is the branch length d.
    // But actually, we don't need to worry about that either because we are taking the proposal density equal to the
    // prior as described below.
    bl_proposer->propose_branches(pp, rng);

    // We want to have:
    // w_r(s_r) = \frac{\gamma^*_r(s_r)}{\gamma^*_{r-1}(s_{r-1})} \frac{1}{\nu^+(s_{r-1} \rightarrow s_r)} (*).
    // Note that \gamma^* is equal to the prior p times the likelihood l.
    // If we take the proposal density equal to the prior, as we have assumed, then the proposal density will be
    // proportional to the ratio of the prior densities. That is,
    // \nu^+(s_{r-1} \rightarrow s_r) = \frac{p(s_r)}{p(s_{r-1})}
    // Thus (*) has some cancellation, becoming
    // w_r(s_r) = \frac{l(s_r)}{l(s_{r-1})};
    // the log of wheich we have here.
    pp->log_likelihood = log_likelihood->calculate_log_likelihood(pp);
    p_from.SetLogWeight(pp->log_likelihood - prev_ll);
    pp->forward_log_density = std::log(fwd_density);
    pp->backward_log_density = std::log(back_density);

    p_from.AddToLogWeight(pp->forward_log_density);
    p_from.AddToLogWeight(-pp->backward_log_density);

    return 0;
}
} // namespace moves
} // namespace sts

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
    auto calc = log_likelihood.get_calculator();
    particle::Particle *part = p_from.GetValuePointer();
    particle::Particle pp = std::make_shared<particle::State>();
    pp->predecessor = *part;
    *part = pp;
    std::vector<particle::Node_ptr> prop_vector = util::uncoalesced_nodes(pp, log_likelihood.get_leaves());

    double prev_ll = log_likelihood(*part);

    // ** First step: perform a uniformly selected merge.
    // Pick two nodes from the prop_vector to join.
    int n1 = rng->UniformDiscrete(0, prop_vector.size() - 1);
    int n2 = rng->UniformDiscrete(0, prop_vector.size() - 2);
    // The following gives the uniform distribution on legal choices that are not n1. Think of taking the uniform
    // distribution on [0,n-2], breaking it at n1 and moving the right hand bit one to the right.
    if(n2 >= n1) n2++;
    pp->node = std::make_shared<particle::Node>(calc);

    // Draw branch lengths.
    pp->node->child1 = std::make_shared<particle::Edge>(prop_vector[n1]);
    pp->node->child2 = std::make_shared<particle::Edge>(prop_vector[n2]);

    // Because the proposal distribution is uniform on merges, the only part of the proposal distribution we have to
    // worry about is the branch length d.
    // But actually, we don't need to worry about that either because we are taking the proposal density equal to the
    // prior as described below.
    bl_proposal(*part, rng);

    // We want to have:
    // w_r(s_r) = \frac{\gamma^*_r(s_r)}{\gamma^*_{r-1}(s_{r-1})} \frac{1}{\nu^+(s_{r-1} \rightarrow s_r)} (*).
    // Note that \gamma^* is equal to the prior p times the likelihood l.
    // If we take the proposal density equal to the prior, as we have assumed, then the proposal density will be
    // proportional to the ratio of the prior densities. That is,
    // \nu^+(s_{r-1} \rightarrow s_r) = \frac{p(s_r)}{p(s_{r-1})}
    // Thus (*) has some cancellation, becoming
    // w_r(s_r) = \frac{l(s_r)}{l(s_{r-1})};
    // the log of wheich we have here.
    (*part)->log_likelihood = log_likelihood(*part);
    p_from.SetLogWeight((*part)->log_likelihood - prev_ll);
    (*part)->forward_log_density = 0.0;

    // Next we multiply by \f$ \nu^-(s_r \rightarrow s_{r-1}) \f$ so that we can correct for multiplicity of particle
    // observation. We can think of this as being the inverse of the number of ways we can get to the current particle
    // s_r from the previous sample. We can think of this as the number of ways to disassemble s_r into something with
    // rank r-1, which is the number of trees in the forest, omitting trees consisting of a single leaf.

    // Count trees (naked leaves don't count) in forest *before* this merge.
    int tc = util::uncoalesced_count_trees(prop_vector);

    // Correct tc for the new number of trees in the forest after this merge.
    if(pp->node->child1->node->is_leaf() && pp->node->child2->node->is_leaf()) {
        // Merged two leaves; number of trees in forest increases by one.
        tc++;
    } else if(!pp->node->child1->node->is_leaf() && !pp->node->child2->node->is_leaf()) {
        // Merged trees; number of trees in forest decreases by one.
        tc--;
    } else {
        // Merged tree + leaf; number of trees unchanged.
        assert(pp->node->child1->node->is_leaf() != pp->node->child2->node->is_leaf());
    }

    (*part)->backward_log_density = -std::log(tc);
    if(tc > 1)
        p_from.AddToLogWeight((*part)->backward_log_density);
    return 0;
}
} // namespace moves
} // namespace sts

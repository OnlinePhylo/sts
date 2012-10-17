#include "phylomoves.hh"
#include <cassert>
#include <vector>
#include <iostream>
#include <cmath>

// forest_likelihood
///The function corresponding to the log likelihood of a forest at specified time and position (up to normalisation)

///  \param X    The state to consider
double forest_likelihood::operator()(const particle& X) const
{
    // Walk backwards through the forest to calculate likelihoods of each tree.
    std::vector<bool> visited;
    double ll_sum = 0;
    std::shared_ptr< phylo_particle > cur = X.pp;
    while(cur != NULL && cur->node != NULL) {
        if(visited.size() < cur->node->id || !visited[ cur->node->id ]) {
            ll_sum += calc->calculate_ll(cur->node, visited);
        }
        cur = cur->predecessor;
    }
    // add background freqs for all uncoalesced leaves
    for(int i = 0; i < leaf_nodes.size(); i++) {
        if(visited.size() > i && visited[i]) continue;
        ll_sum += calc->calculate_ll(leaf_nodes[i], visited);
    }

    return ll_sum;
}

/// Get the calculator
inline std::shared_ptr<online_calculator> forest_likelihood::get_calculator() const { return calc; }

/// Get the vector representing \\perp
inline const std::vector<std::shared_ptr<phylo_node>> forest_likelihood::get_leaves() const { return leaf_nodes; }
// /forest_likelihood

// mcmc_move
/// Function call for use with smctc - calls user-defined do_move, tracks result.

/// \param time Generation number
/// \param from Source particle
/// \param rng Random number source
int mcmc_move::operator()(long time, smc::particle<particle>& from, smc::rng *rng)
{
    attempted++;
    int result = do_move(time, from, rng);
    if(result) accepted++;
    return result;
}
// /mcmc_move

// uniform_bl_mcmc_move
/// Uniform change to branch lengths for the current node

/// Change the branch lengths for the current node by drawing from a uniform distribution between -amount and amount.
///  \param time  generation number
///  \param from  Source particle
///  \param rng   Random number source
int uniform_bl_mcmc_move::do_move(long time, smc::particle<particle>& from, smc::rng* rng) const
{
    std::shared_ptr<online_calculator> calc = log_likelihood.get_calculator();
    particle* part = from.GetValuePointer();
    std::shared_ptr< phylo_node > cur_node = part->pp->node;
    std::shared_ptr< phylo_node > new_node = std::make_shared<phylo_node>(calc);
    new_node->child1 = cur_node->child1;
    new_node->child2 = cur_node->child2;
    new_node->id = calc->get_id();

    double cur_ll = log_likelihood(*part);
    // Choose an amount to shift the node height uniformly at random.
    double shift = rng->Uniform(-amount, amount);
    // If the shift amount would create a negative node height we will reflect it back to a positive number.
    // This means the proposal density for the reflection zone is double the rest of the area, but the back-proposal
    // probability is also double in the same area so these terms cancel in the Metropolis-Hastings ratio.

    // Now calculate the new node heights - shift both heights for now: ultrametric
    new_node->child1->length = fabs(new_node->child1->length + shift);
    new_node->child2->length = fabs(new_node->child2->length + shift);
    part->pp->node = new_node;

    double alpha = exp(log_likelihood(*part) - cur_ll);
    if(alpha < 1 && rng->UniformS() > alpha) {
        // Move rejected, restore the original node.
        part->pp->node = cur_node;
        return false;
    }
    // Accept the new state.
    return true;
}
// /uniform_bl_mcmc_move

// smc_move
int smc_move::operator()(long t, smc::particle<particle>& p, smc::rng* r)
{
    call_count++;
    return do_move(t, p, r);
}
// /smc_move

// fmove

/// Perform a particle move.
/// It merges trees in p_from in place, and updates their weights.
/// \f[
/// w_r(s_r) = \frac{\gamma_r(s_r)}{\gamma_{r-1}(s_{r-1})} \frac{\nu^-(s_r \rightarrow s_{r-1})}{\nu^+(s_{r-1} \rightarrow s_r)}.
/// \f]

///\param time The sampler iteration.
///\param p_from The particle to move.
///\param rng  A random number generator.
int rooted_merge::do_move(long time, smc::particle<particle>& p_from, smc::rng* rng) const
{
    std::shared_ptr<online_calculator> calc = log_likelihood.get_calculator();
    particle* part = p_from.GetValuePointer();
    std::shared_ptr< phylo_particle > pp = std::make_shared<phylo_particle>();
    pp->predecessor = part->pp;
    part->pp = pp;
    std::vector< std::shared_ptr< phylo_node > > prop_vector = uncoalesced_nodes(pp, log_likelihood.get_leaves());

    double prev_ll = log_likelihood(*part);

    // ** First step: perform a uniformly selected merge.
    // Pick two nodes from the prop_vector to join.
    int n1 = rng->UniformDiscrete(0, prop_vector.size() - 1);
    int n2 = rng->UniformDiscrete(0, prop_vector.size() - 2);;
    if(n2 >= n1) n2++;
    pp->node = std::make_shared<phylo_node>(calc);
    pp->node->id = calc->get_id();

    // Draw branch lengths.
    // For ultrametric case, need to set d1 = d2
    // double d1 = pRng->Exponential(1.0), d2 = pRng->Exponential(1.0);
    // NOTE: if the mean of this distribution is changed from 1.0 then we will need to also update the formula for
    // d_prob below.
    double d1 = rng->Exponential(1.0), d2 = d1;
    pp->node->child1 = std::make_shared<edge>(prop_vector[n1], d1);
    pp->node->child2 = std::make_shared<edge>(prop_vector[n2], d2);
    pp->node->calc_height();

    // ** Second step: calculate weights.
    // d_prob is q(s->s'), *not* on log scale
    // Ultrametric for now.
    // double d_prob = exp(-d1 - d2);
    // Because the proposal distribution is uniform on merges, the only part of the proposal distribution we have to
    // worry about is the branch length d. Because we are drawing it from an exponential distribution we have:
    double d_prob = exp(-d1);

    // Note: when proposing from exponential(1.0) the below can be simplified to just adding d1 and d2
    // We want to have:
    // w_r(s_r) = \frac{\gamma_r(s_r)}{\gamma_{r-1}(s_{r-1})} \frac{1}{\nu^+(s_{r-1} \rightarrow s_r)}.
    p_from.SetLogWeight(log_likelihood(*part) - prev_ll - log(d_prob));

    // Next we multiply by \f$ \nu^-(s_r \rightarrow s_{r-1}) \f$ so that we can correct for multiplicity of particle
    // observation. We can think of this as being the inverse of the number of ways we can get to the current particle
    // s_r from the previous sample. We can think of this as the number of ways to disassemble s_r into something with
    // rank r-1, which is the number of trees in the forest, omitting trees consisting of a single leaf.
    // We add one below because prop_vector is from the previous particle, which had one less tree than the post-merging
    // particle does.
    const int tc = tree_count(prop_vector)+1;
    if(tc > 1)
        p_from.AddToLogWeight(-log(tc));
    return 0;
}
// /fmove

// smc_init
///A function to initialise particles

/// \param rng A pointer to the random number generator which is to be used
smc::particle<particle> smc_init::operator()(smc::rng* rng)
{
    particle value;
    // initial particles have all sequences uncoalesced
    value.pp = std::make_shared< phylo_particle >();
    // Note that the likelihood of the equivalent \perp particles doesn't matter. We set it to zero.
    return smc::particle<particle>(value, 0.);
}
// /smc_init

#include "phylomoves.hh"
#include <vector>

// forest_likelihood
///The function corresponding to the log likelihood of a forest at specified time and position (up to normalisation)

///  \param calc The likelihood calculator
///  \param time The current time (i.e. the number of coalescence events so far)
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

inline std::shared_ptr<OnlineCalculator> forest_likelihood::get_calculator() const { return calc; }
inline const std::vector<std::shared_ptr<phylo_node>> forest_likelihood::get_leaves() const { return leaf_nodes; }
// /forest_likelihood

// mcmc_move
/// Function call for use with smctc - calls user-defined do_move, tracks result.
int mcmc_move::operator()(long time, smc::particle<particle>& from, smc::rng *rng) {
    attempted++;
    int result = do_move(time, from, rng);
    if(result) accepted++;
    return result;
}
// /mcmc_move

// uniform_bl_mcmc_move
/// Change the branch lengths for the current node by drawing from a uniform distribution between -amount and amount.
int uniform_bl_mcmc_move::do_move(long time, smc::particle<particle>& from, smc::rng* rng) const {
    std::shared_ptr<OnlineCalculator> calc = log_likelihood.get_calculator();
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
    new_node->child1->length = abs(new_node->child1->length + shift);
    new_node->child2->length = abs(new_node->child2->length + shift);
    new_node->calc_height();
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
int smc_move::operator()(long t, smc::particle<particle>& p, smc::rng* r) {
    call_count++;
    return do_move(t, p, r);
}


// /smc_move

// fmove

///The proposal function.

///\param time The sampler iteration.
///\param p_from The particle to move.
///\param rng  A random number generator.
int rooted_merge::do_move(long time, smc::particle<particle>& p_from, smc::rng* rng) const
{
    std::shared_ptr<OnlineCalculator> calc = log_likelihood.get_calculator();
    particle* part = p_from.GetValuePointer();
    std::shared_ptr< phylo_particle > pp = std::make_shared<phylo_particle>();
    pp->predecessor = part->pp;
    part->pp = pp;
    std::vector< std::shared_ptr< phylo_node > > prop_vector = uncoalesced_nodes(pp, log_likelihood.get_leaves());

    // Pick two nodes from the prop_vector to join.
    int n1 = rng->UniformDiscrete(0, prop_vector.size() - 1);
    int n2 = rng->UniformDiscrete(0, prop_vector.size() - 2);;
    if(n2 >= n1) n2++;
    pp->node = std::make_shared<phylo_node>(calc);
    //pp->node = std::make_shared< phylo_node >();
    pp->node->id = calc->get_id();

    // Propose branch lengths.
    // For ultrametric case, need to set d1 = d2
    // double d1 = pRng->Exponential(1.0), d2 = pRng->Exponential(1.0);
    double d1 = rng->Exponential(1.0), d2 = d1;
    pp->node->child1 = std::make_shared<edge>(prop_vector[n1], d1);
    pp->node->child2 = std::make_shared<edge>(prop_vector[n2], d2);
    pp->node->calc_height();

    // d_prob is q(s->s'), *not* on log scale
    // Ultrametric for now
    //double d_prob = exp(-d1 - d2);
    double d_prob = exp(-d1);

    // Note: when proposing from exponential(1.0) the below can be simplified to just adding d1 and d2
    p_from.AddToLogWeight(log_likelihood(*part) - log(d_prob));

    // Add reverse transition probability q(r' -> r)
    // 1/(# of trees in forest), omitting trees consisting of a single leaf
    // This prevents over-counting
    const int tc = tree_count(prop_vector);
    if(tc > 1)
        p_from.AddToLogWeight(-log(tc));
    return 0;
}
// /fmove

// smc_init
///A function to initialise particles

/// \param rng A pointer to the random number generator which is to be used
smc::particle<particle> smc_init::operator()(smc::rng* rng) {
    particle value;
    // initial particles have all sequences uncoalesced
    value.pp = std::make_shared< phylo_particle >();
    // loglike should just be the background distribution on character state frequencies
    return smc::particle<particle>(value, log_likelihood(value));
}
// /smc_init

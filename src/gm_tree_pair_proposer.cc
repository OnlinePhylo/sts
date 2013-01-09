#include "gm_tree_pair_proposer.h"
#include "util.h"
#include "node.h"

#include <cassert>
#include <iterator>
#include <vector>
#include <iostream>
#include <utility>

using namespace std;
using sts::particle::Node_ptr;

namespace sts
{
namespace moves
{

GM_tree_pair_proposer::GM_tree_pair_proposer(sts::likelihood::Forest_likelihood* log_likelihood, const sts::guidedmerge::GM_tree tree, double poisson_mean) :
    log_likelihood(log_likelihood),
    mu(poisson_mean)
{
    gm_trees[nullptr] = tree;
}

void GM_tree_pair_proposer::operator()(particle::Particle pp, smc::rng* rng, particle::Node_ptr& a, particle::Node_ptr& b, double& fwd_density, double& back_density)
{

    // propose a join uniformly at random
    std::vector<particle::Node_ptr> prop_vector = util::uncoalesced_nodes(pp, log_likelihood->get_leaves());

    //fwd_density = 1.0 / ((prop_vector.size() - 1) * (prop_vector.size() - 2));
    // Since merges are uniform, setting to 1
    fwd_density = 1.0;
    back_density = 1.0;

    // Nothing merged yet
    if(prop_vector.size() == 2) { // last merge
        a = prop_vector[0];
        b = prop_vector[1];
        return;
    }

    size_t n_tries = 0;
    // TODO: something smarter
    while(++n_tries < 1000) {
        unsigned int k = rng->Poisson(mu);

        // pp is the current node being worked on
        // In the first generation, its predecessor is null, for which the GM_tree is assigned in the constructor
        auto merges = gm_trees.at(pp->predecessor->node.get()).find_k_distance_merges(k);
        if(!merges.empty()) {
            vector<pair<Node_ptr,Node_ptr>> v(begin(merges), end(merges));
            size_t i = rng->UniformDiscrete(0, v.size() - 1);
            a = v[i].first;
            b = v[i].second;
            fwd_density = gsl_ran_poisson_pdf(k, mu);

            // Update gm_tree
            gm_trees[pp->node.get()] = gm_trees[pp->predecessor->node.get()];
            assert(pp->node != nullptr);
            gm_trees[pp->node.get()].merge(a, b, pp->node);
            break;
        }
    }

}


} // namespace moves
} // namespace sts

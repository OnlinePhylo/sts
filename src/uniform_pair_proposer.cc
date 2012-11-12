#include <vector>
#include "uniform_pair_proposer.h"
#include "util.h"

using namespace std;

namespace sts
{
namespace moves
{

void Uniform_pair_proposer::operator()(particle::Particle pp, smc::rng* rng, particle::Node_ptr& a, particle::Node_ptr& b, double& fwd_density, double& back_density)
{
  // propose a join uniformly at random
    auto calc = log_likelihood.get_calculator();
    std::vector<particle::Node_ptr> prop_vector = util::uncoalesced_nodes(pp, log_likelihood.get_leaves());

    // ** First step: perform a uniformly selected merge.
    // Pick two nodes from the prop_vector to join.
    int n1 = rng->UniformDiscrete(0, prop_vector.size() - 1);
    int n2 = rng->UniformDiscrete(0, prop_vector.size() - 2);
    // The following gives the uniform distribution on legal choices that are not n1. Think of taking the uniform
    // distribution on [0,n-2], breaking it at n1 and moving the right hand bit one to the right.
    if(n2 >= n1) n2++;
    a = prop_vector[n1];
    b = prop_vector[n2];
    
    fwd_density = 1.0 / ((prop_vector.size() - 1) * (prop_vector.size() - 2));

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
    
    back_density = 1.0;
    if(tc > 1)
      back_density = tc;
}


} // namespace moves
} // namespace sts

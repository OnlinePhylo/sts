#include "child_swap_mcmc_move.h"
#include "edge.h"
#include "util.h"

#include <cassert>
#include <cmath>

namespace sts
{
namespace moves
{

int Child_swap_mcmc_move::do_move(long time, smc::particle<particle::Particle>& from, smc::rng* rng) const
{
    // Don't make a move for leaf nodes.
    if((*from.GetValuePointer())->node->is_leaf()) return 0;

    auto calc = log_likelihood.get_calculator();
    particle::Particle part = *from.GetValuePointer();
    particle::Node_ptr cur_node = part->node;
    particle::Node_ptr new_node = std::make_shared<particle::Node>(*cur_node);

    double cur_ll = log_likelihood(part);

    std::vector<particle::Node_ptr> prop_vector = util::uncoalesced_nodes(pp, log_likelihood.get_leaves());

    // Choose one of the two children randomly
    const int child_idx = rng->UniformDiscrete(0, 1);

    std::shared_ptr<particle::Edge> *e;

    if(child_idx == 0) {
        e = &(new_node->child1);
    } else if(child_idx == 1) {
        e = &(new_node->child2);
    } else {
        assert(false);
    }

    // Choose an uncoalesced node to swap, retaining current length
    const int n = rng->UniformDiscrete(0, prop_vector.size() - 1);

    // TODO: Should we do something to set edge length?
    e->reset(new particle::Edge(prop_vector[n], (*e)->length));

    part->node = new_node;

    double alpha = exp(log_likelihood(part) - cur_ll);
    if(alpha < 1 && rng->UniformS() > alpha) {
        // Move rejected, restore the original node.
        part->node = cur_node;
        new_node.reset();
        return false;
    }
    // Accept the new state.
    return true;
}

}
}
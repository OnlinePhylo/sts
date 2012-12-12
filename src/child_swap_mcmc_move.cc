#include "child_swap_mcmc_move.h"
#include "edge.h"
#include "util.h"

#include <cassert>
#include <cmath>

namespace sts
{
namespace moves
{

/// Propose a swap between one of the two children of this node and an uncoalesced node.
void Child_swap_mcmc_move::propose_move(long time, particle::Particle& part, smc::rng* rng) const
{
    // Ignore leaf nodes
    if(part->node->is_leaf()) return;

    std::vector<particle::Node_ptr> prop_vector = util::uncoalesced_nodes(part, log_likelihood.get_leaves());

    // Choose one of the two children randomly
    const int child_idx = rng->UniformDiscrete(0, 1);

    std::shared_ptr<particle::Edge> *e;

    if(child_idx == 0) {
        e = &(part->node->child1);
    } else if(child_idx == 1) {
        e = &(part->node->child2);
    } else {
        assert(false);
    }

    // Choose an uncoalesced node to swap, retaining current length
    const int n = rng->UniformDiscrete(0, prop_vector.size() - 1);

    // TODO: Should we do something to set edge length?
    e->reset(new particle::Edge(prop_vector[n], (*e)->length));
}

}
}

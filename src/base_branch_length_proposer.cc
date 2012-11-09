#include "base_branch_length_proposer.h"
#include "node.h"
#include "state.h"
#include "edge.h"
#include <cassert>

namespace sts
{
namespace moves
{

double Base_branch_length_proposer::operator()(particle::Particle part, smc::rng *rng)
{
    Branch_lengths p = propose(part, rng); // This is where the subclassing action happens.
    std::shared_ptr<particle::Node> node = part->node;

    // Children should be initialized
    assert(node->child1);
    assert(node->child2);
    node->child1->length = p.first;
    node->child2->length = p.second;
    return log_proposal_density(p.first) + log_proposal_density(p.second);
}

Base_branch_length_proposer::Branch_lengths Base_branch_length_proposer::propose(particle::Particle part, smc::rng *rng)
{
    double d1 = propose_bl(rng), d2 = propose_bl(rng);
    return Base_branch_length_proposer::Branch_lengths(d1, d2);
}

} // namespace moves
} // namespace sts

#ifndef STS_MOVES_ROOTED_MERGE_HPP
#define STS_MOVES_ROOTED_MERGE_HPP
#include <cassert>
#include <cmath>
#include <iostream>
#include <functional>
#include <memory>
#include <vector>

#include "sts/likelihood/forest_likelihood.hpp"
#include "sts/moves/exponential_branch_length_proposer.hpp"
#include "sts/moves/smc_move.hpp"
#include "sts/particle/phylo_node.hpp"
#include "sts/particle/phylo_particle.hpp"
#include "sts/util.hpp"

namespace sts
{
namespace moves
{

/// \class rooted_merge
/// \brief Merge of two nodes, with exponential branch length proposal
class rooted_merge: public smc_move
{
public:
    /// Branch length proposal function.
    /// Accepts two parameters: a phylo_node with initialized edges and a random source; returns the log-likelihood.
    typedef std::function<double(particle::particle, smc::rng*)> bl_proposal_fn;

    /// Constructor

    /// Initializes with exponential_branch_length_proposal with mean 1.0.
    explicit rooted_merge(sts::likelihood::forest_likelihood& log_likelihood) : smc_move(log_likelihood),
        bl_proposal(exponential_branch_length_proposer(1.0)) {};

    /// Constructor

    /// \param bl_proposal Source of branch length proposals.
    rooted_merge(sts::likelihood::forest_likelihood& log_likelihood,
                 bl_proposal_fn bl_proposal) : smc_move(log_likelihood), bl_proposal(bl_proposal) {};

    int do_move(long, smc::particle<particle::particle>&, smc::rng*) const;

protected:
    /// Branch length proposal generator
    bl_proposal_fn bl_proposal;
};

// Implementation

/// Perform a particle move.

/// It merges trees in \c p_from in place, and updates their weights.
/// \f[
/// w_r(s_r) = \frac{\gamma_r(s_r)}{\gamma_{r-1}(s_{r-1})} \frac{\nu^-(s_r \rightarrow s_{r-1})}{\nu^+(s_{r-1} \rightarrow s_r)}.
/// \f]

///\param time The sampler iteration.
///\param p_from The particle to move.
///\param rng  A random number generator.
int rooted_merge::do_move(long time, smc::particle<particle::particle>& p_from, smc::rng* rng) const
{
    auto calc = log_likelihood.get_calculator();
    particle::particle *part = p_from.GetValuePointer();
    particle::particle pp = std::make_shared<particle::phylo_particle>();
    pp->predecessor = *part;
    *part = pp;
    std::vector<particle::node> prop_vector = util::uncoalesced_nodes(pp, log_likelihood.get_leaves());

    double prev_ll = log_likelihood(*part);

    // ** First step: perform a uniformly selected merge.
    // Pick two nodes from the prop_vector to join.
    int n1 = rng->UniformDiscrete(0, prop_vector.size() - 1);
    int n2 = rng->UniformDiscrete(0, prop_vector.size() - 2);;
    if(n2 >= n1) n2++;
    pp->node = std::make_shared<particle::phylo_node>(calc);

    // Draw branch lengths.
    // For ultrametric case, need to set d1 = d2
    // double d1 = pRng->Exponential(1.0), d2 = pRng->Exponential(1.0);
    // NOTE: if the mean of this distribution is changed from 1.0 then we will need to also update the formula for
    // d_prob below.
    pp->node->child1 = std::make_shared<particle::edge>(prop_vector[n1]);
    pp->node->child2 = std::make_shared<particle::edge>(prop_vector[n2]);

    // Because the proposal distribution is uniform on merges, the only part of the proposal distribution we have to
    // worry about is the branch length d.
    // This is returned from the proposal function.
    double bl_ll = bl_proposal(*part, rng);
    pp->node->calc_height();

    // We want to have:
    // w_r(s_r) = \frac{\gamma_r(s_r)}{\gamma_{r-1}(s_{r-1})} \frac{1}{\nu^+(s_{r-1} \rightarrow s_r)}.
    p_from.SetLogWeight(log_likelihood(*part) - prev_ll - bl_ll);

    // Next we multiply by \f$ \nu^-(s_r \rightarrow s_{r-1}) \f$ so that we can correct for multiplicity of particle
    // observation. We can think of this as being the inverse of the number of ways we can get to the current particle
    // s_r from the previous sample. We can think of this as the number of ways to disassemble s_r into something with
    // rank r-1, which is the number of trees in the forest, omitting trees consisting of a single leaf.
    int tc = util::tree_count(prop_vector);
    if(pp->node->child1->node->is_leaf() == pp->node->child2->node->is_leaf()) {
        // When the merge is between two leaf nodes, or two non-leaf nodes, this merge increases the tree count by 1.
        // A merge between a leaf and a non-leaf node does not change the tree count.
        tc++;
    }

    if(tc > 1)
        p_from.AddToLogWeight(-log(tc));
    return 0;
}
} // namespace moves
} // namespace sts

#endif // STS_MOVES_ROOTED_MERGE_HPP

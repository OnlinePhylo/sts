#ifndef STS_GM_TREE_PAIR_PROPOSER_H
#define STS_GM_TREE_PAIR_PROPOSER_H

#include "forest_likelihood.h"
#include "gm_tree.h"
#include "node_ptr.h"
#include "state.h"

#include <smctc.hh>
#include <string>
#include <vector>
#include <unordered_map>

namespace sts
{
namespace moves
{

/// \class GM_tree_pair_proposer
/// \brief "Propose" pairs of nodes against a GM tree
class GM_tree_pair_proposer
{
public:
    /// Instantiate
    GM_tree_pair_proposer(sts::likelihood::Forest_likelihood* log_likelihood, const sts::guidedmerge::GM_tree tree, double poisson_mean);
    void operator()(particle::Particle, smc::rng*, particle::Node_ptr& a, particle::Node_ptr& b, double& fwd_density, double& back_density);
private:
    sts::likelihood::Forest_likelihood* log_likelihood;
    std::unordered_map<const sts::particle::Node*,sts::guidedmerge::GM_tree> gm_trees;
    /// \brief Mean of the poisson distribution
    double mu;
};


} // namespace moves
} // namespace sts

#endif // STS_GM_TREE_PAIR_PROPOSER_H

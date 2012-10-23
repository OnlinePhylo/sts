#ifndef STS_LIKELIHOOD_FOREST_LIKELIHOOD_HPP
#define STS_LIKELIHOOD_FOREST_LIKELIHOOD_HPP

#include <memory>
#include <vector>

#include "sts/particle/phylo_node.hpp"
#include "sts/particle/phylo_particle.hpp"
#include "sts/likelihood/online_calculator.hpp"

namespace sts
{
namespace likelihood
{

/// Class to calculate the likelihood of a forest
class forest_likelihood
{
public:
    /// Constructor

    ///  \param calc Initialized likelihood calculator
    ///  \param leaf_nodes Vector representing \\perp
    explicit forest_likelihood(std::shared_ptr<online_calculator> calc,
                               std::vector<std::shared_ptr<particle::phylo_node>> leaf_nodes) : calc(calc), leaf_nodes(leaf_nodes) {};
    /// Copy constructor
    explicit forest_likelihood(const forest_likelihood &other) : calc(other.calc), leaf_nodes(other.leaf_nodes) {};

    double operator()(const particle::particle&) const;

    const std::vector<std::shared_ptr<particle::phylo_node>> get_leaves() const;
    std::shared_ptr<online_calculator> get_calculator() const;
private:
    std::shared_ptr<online_calculator> calc;
    std::vector<std::shared_ptr<particle::phylo_node>> leaf_nodes;
};

///The function corresponding to the log likelihood of a forest at specified time and position (up to normalisation)

///  \param X    The state to consider
double forest_likelihood::operator()(const particle::particle& X) const
{
    // Walk backwards through the forest to calculate likelihoods of each tree.
    std::vector<bool> visited;
    double ll_sum = 0;
    particle::particle cur = X;
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
inline const std::vector<std::shared_ptr<particle::phylo_node>> forest_likelihood::get_leaves() const { return leaf_nodes; }

} // namespace likelihood
} // namespace sts

#endif // STS_LIKELIHOOD_FOREST_LIKELIHOOD_HPP

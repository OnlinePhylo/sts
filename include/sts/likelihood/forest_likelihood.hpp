#ifndef STS_LIKELIHOOD_FOREST_LIKELIHOOD_HPP
#define STS_LIKELIHOOD_FOREST_LIKELIHOOD_HPP

#include <memory>
#include <vector>

#include "sts/particle/node.hpp"
#include "sts/particle/state.hpp"
#include "sts/likelihood/online_calculator.hpp"

namespace sts
{
namespace likelihood
{

/// Class to calculate the likelihood of a forest
class Forest_likelihood
{
public:
    /// Constructor

    ///  \param calc Initialized likelihood calculator
    ///  \param leaf_nodes Vector representing \\perp
    explicit Forest_likelihood(std::shared_ptr<Online_calculator> calc,
                               std::vector<particle::Node_ptr> leaf_nodes) : calc(calc), leaf_nodes(leaf_nodes) {};
    /// Copy constructor
    explicit Forest_likelihood(const Forest_likelihood &other) : calc(other.calc), leaf_nodes(other.leaf_nodes) {};

    double operator()(const particle::Particle&) const;

    const std::vector<particle::Node_ptr> get_leaves() const;
    std::shared_ptr<Online_calculator> get_calculator() const;
private:
    std::shared_ptr<Online_calculator> calc;
    std::vector<particle::Node_ptr> leaf_nodes;
};

///The function corresponding to the log likelihood of a forest at specified time and position (up to normalisation)

///  \param X    The state to consider
double Forest_likelihood::operator()(const particle::Particle& X) const
{
    // Walk backwards through the forest to calculate likelihoods of each tree.
    std::unordered_set<particle::Node_ptr> visited;
    double ll_sum = 0;
    particle::Particle cur = X;
    while(cur != NULL && cur->node != NULL) {
        if(visited.count(cur->node) == 0) {
            ll_sum += calc->calculate_ll(cur->node, visited);
        }
        cur = cur->predecessor;
    }
    // add background freqs for all uncoalesced leaves
    for(int i = 0; i < leaf_nodes.size(); i++) {
        if(visited.count(leaf_nodes[i]) != 0) continue;
        ll_sum += calc->calculate_ll(leaf_nodes[i], visited);
    }

    return ll_sum;
}

/// Get the calculator
inline std::shared_ptr<Online_calculator> Forest_likelihood::get_calculator() const { return calc; }

/// Get the vector representing \\perp
inline const std::vector<particle::Node_ptr> Forest_likelihood::get_leaves() const { return leaf_nodes; }

} // namespace likelihood
} // namespace sts

#endif // STS_LIKELIHOOD_FOREST_LIKELIHOOD_HPP

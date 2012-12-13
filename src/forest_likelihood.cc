#include "forest_likelihood.h"

#include "online_calculator.h"
#include "node.h"
#include "state.h"

#include <unordered_set>

using namespace std;

namespace sts
{
namespace likelihood
{
/// Shortcut for \c calculate_log_likelihood

///  \param X The state to consider
double Forest_likelihood::operator()(const particle::Particle& X) const
{
    return this->calculate_log_likelihood(X);
}

/// \brief Calculate the log-likelihood of a forest

/// The function corresponding to the log likelihood of a forest at specified time and position (up to normalisation)
/// \param X Particle to act on
double Forest_likelihood::calculate_log_likelihood(const particle::Particle& X) const
{
    // Walk backwards through the forest to calculate likelihoods of each tree.
    unordered_set<particle::Node_ptr> visited;
    double ll_sum = 0;
    particle::Particle cur = X;
    while(cur != nullptr && cur->node != nullptr) {
        if(visited.count(cur->node) == 0) {
            ll_sum += calc->calculate_ll(cur->node, visited);
        }
        cur = cur->predecessor;
    }
    // add background freqs for all uncoalesced leaves
    for(size_t i = 0; i < leaf_nodes.size(); i++) {
        if(visited.count(leaf_nodes[i]) != 0) continue;
        ll_sum += calc->calculate_ll(leaf_nodes[i], visited);
    }

    return ll_sum;
}


/// Get the calculator
shared_ptr<Online_calculator> Forest_likelihood::get_calculator() const { return calc; }

/// Get the vector representing \\perp
const vector<particle::Node_ptr> Forest_likelihood::get_leaves() const { return leaf_nodes; }

} // namespace likelihood
} // namespace sts

#ifndef STS_PARTICLE_PHYLO_PARTICLE_H
#define STS_PARTICLE_PHYLO_PARTICLE_H

#include <memory>
#include <stack>
#include <string>
#include <unordered_map>
#include <Bpp/Phyl/TreeTemplateTools.h>
#include "online_calculator.h"
#include "node_ptr.h"
#include "particle.h"

namespace sts
{
namespace particle
{

/// \class State
/// \brief A forest in the SMC.
/// This class stores the SMC forest implicitly, by specifying the collections
/// of mergers that must be made in order to get the forest from \f$\perp\f$, the
/// completely un-merged state.
class State
{
public:
    State() :
        node(nullptr),
        predecessor(nullptr),
        forward_log_density(0.0),
        backward_log_density(0.0),
        partial_log_likelihood(0.0) {};

    /// The merge novel to this particle. If \c nullptr then the particle is \f$\perp\f$.
    std::shared_ptr<Node> node;
    /// The predecessor particles, which specify the rest of the merges for this particle.
    std::shared_ptr<State> predecessor;

    /// Forward log proposal density: \f$\nu^+(s_{r-1} \rightarrow s_r)\f$
    double forward_log_density;
    /// Backward log proposal density: \f$\nu^-(s_{r} \rightarrow s_{r-1})\f$
    double backward_log_density;
    /// partial likelihood of this state: \f$\frac{L_{s+1}}{L_s}\f$
    // (\f$\frac{\gamma*_{s+1}}{\gamma*_s}\f$).
    double partial_log_likelihood;

    /// Make a State from a bpp Tree
    static std::shared_ptr<State>
    of_tree(std::shared_ptr<likelihood::Online_calculator>, bpp::TreeTemplate<bpp::Node> &, std::unordered_map<sts::particle::Node_ptr, std::string>&);

    /// Make a State from a Newick tree string
    static std::shared_ptr<State>
    of_newick_string(std::shared_ptr<likelihood::Online_calculator>, std::string &, std::unordered_map<sts::particle::Node_ptr, std::string>&);
};

} // namespace particle
} // namespace sts

#endif // STS_PARTICLE_PHYLO_PARTICLE_H
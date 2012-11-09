#ifndef STS_PARTICLE_PHYLO_PARTICLE_HPP
#define STS_PARTICLE_PHYLO_PARTICLE_HPP

#include <memory>
#include <stack>
#include <string>
#include <unordered_map>
#include <Bpp/Phyl/TreeTemplateTools.h>
#include "online_calculator.h"

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
    /// The merge novel to this particle. If NULL then the particle is \f$\perp\f$.
    std::shared_ptr<Node> node;
    /// The predecessor particles, which specify the rest of the merges for this particle.
    std::shared_ptr<State> predecessor;

    /// Make a State from a bpp Tree
    static std::shared_ptr<State>
    of_tree(std::shared_ptr<likelihood::Online_calculator>, bpp::TreeTemplate<bpp::Node> &, std::unordered_map<sts::particle::Node_ptr, std::string>&);

    /// Make a State from a Newick tree string
    static std::shared_ptr<State>
    of_newick_string(std::shared_ptr<likelihood::Online_calculator>, std::string &, std::unordered_map<sts::particle::Node_ptr, std::string>&);
};

} // namespace particle
} // namespace sts

#endif // STS_PARTICLE_PHYLO_PARTICLE_HPP


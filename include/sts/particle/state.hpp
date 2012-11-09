#ifndef STS_PARTICLE_PHYLO_PARTICLE_HPP
#define STS_PARTICLE_PHYLO_PARTICLE_HPP

#include <memory>
#include <stack>
#include <string>
#include <unordered_map>
#include <Bpp/Phyl/TreeTemplateTools.h>
#include "sts/likelihood/online_calculator.hpp"

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
    of_tree(std::shared_ptr<likelihood::online_calculator>, bpp::TreeTemplate<bpp::Node> &, std::unordered_map<sts::particle::node_ptr, std::string>&);

    /// Make a State from a Newick tree string
    static std::shared_ptr<State>
    of_newick_string(std::shared_ptr<likelihood::online_calculator>, std::string &, std::unordered_map<sts::particle::node_ptr, std::string>&);
};


/// A particle in the SMC
typedef std::shared_ptr<State> particle;

particle State::of_tree(std::shared_ptr<likelihood::online_calculator> calc, bpp::TreeTemplate<bpp::Node> &tree,
                                 std::unordered_map<sts::particle::node_ptr, std::string>& names)
{
    particle p = std::make_shared<State>();
    p->node = Node::of_tree(calc, tree, names);
    if(p->node->is_leaf())
        return p;

    particle prev = p;
    std::stack<sts::particle::node_ptr> node_stack;
    node_stack.push(p->node->child1->node);
    node_stack.push(p->node->child2->node);
    while(!node_stack.empty()) {
        particle cur = std::make_shared<State>();
        cur->node = node_stack.top();
        node_stack.pop();
        prev->predecessor = cur;
        if(!cur->node->is_leaf()) {
            node_stack.push(cur->node->child1->node);
            node_stack.push(cur->node->child2->node);
        }
        prev = cur;
    }

    return p;
}


particle State::of_newick_string(std::shared_ptr<likelihood::online_calculator> calc, std::string &tree_string,
        std::unordered_map<sts::particle::node_ptr, std::string>& names)
{
    bpp::TreeTemplate<bpp::Node> *tree = bpp::TreeTemplateTools::parenthesisToTree(tree_string);
    particle node = State::of_tree(calc, *tree, names);
    delete tree;
    return node;
}

} // namespace particle
} // namespace sts

#endif // STS_PARTICLE_PHYLO_PARTICLE_HPP


#ifndef STS_PARTICLE_PHYLO_PARTICLE_HPP
#define STS_PARTICLE_PHYLO_PARTICLE_HPP

#include <memory>
#include <stack>
#include <Bpp/Phyl/TreeTemplateTools.h>
#include "sts/likelihood/online_calculator.hpp"

namespace sts
{
namespace particle
{
/// \class phylo_particle
/// A forest in the SMC.
///
/// This class stores the SMC forest implicitly, by specifying the collections
/// of mergers that must be made in order to get the forest from \perp, the
/// completely un-merged state.
class phylo_particle
{
public:
    // The merge novel to this particle. If NULL then the particle is \perp.
    std::shared_ptr<phylo_node> node;
    // The predecessor particles, which specify the rest of the merges for this particle.
    std::shared_ptr<phylo_particle> predecessor;

    /// Make a phylo_particle from a bpp Tree
    static std::shared_ptr<phylo_particle>
    of_tree(std::shared_ptr<likelihood::online_calculator>, bpp::TreeTemplate<bpp::Node> &);

    /// Make a phylo_particle from a Newick tree string
    static std::shared_ptr<phylo_particle>
    of_newick_string(std::shared_ptr<likelihood::online_calculator>, std::string &);
};

class particle
{
public:
    std::shared_ptr<phylo_particle> pp;
};

std::shared_ptr<phylo_particle> phylo_particle::of_tree(std::shared_ptr<likelihood::online_calculator> calc, bpp::TreeTemplate<bpp::Node> &tree)
{
    std::shared_ptr<phylo_particle> particle = std::make_shared<phylo_particle>();
    particle->node = phylo_node::of_tree(calc, tree);
    if(particle->node->is_leaf())
        return particle;

    std::shared_ptr<phylo_particle> prev = particle;
    std::stack<std::shared_ptr<phylo_node>> node_stack;
    node_stack.push(particle->node->child1->node);
    node_stack.push(particle->node->child2->node);
    while(!node_stack.empty()) {
        std::shared_ptr<phylo_particle> cur = std::make_shared<phylo_particle>();
        cur->node = node_stack.top();
        node_stack.pop();
        prev->predecessor = cur;
        if(!cur->node->is_leaf()) {
            node_stack.push(cur->node->child1->node);
            node_stack.push(cur->node->child2->node);
        }
        prev = cur;
    }

    return particle;
}


std::shared_ptr <phylo_particle> phylo_particle::of_newick_string(std::shared_ptr<likelihood::online_calculator> calc, std::string &tree_string)
{
    bpp::TreeTemplate<bpp::Node> *tree = bpp::TreeTemplateTools::parenthesisToTree(tree_string);
    std::shared_ptr<phylo_particle> node = phylo_particle::of_tree(calc, *tree);
    delete tree;
    return node;
}

} // namespace particle
} // namespace sts

#endif // STS_PARTICLE_PHYLO_PARTICLE_HPP

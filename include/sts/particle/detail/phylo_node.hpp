/// \file detail/phylo_node.hpp

/// \brief phylo_node definition.

#ifndef STS_PARTICLE_DETAIL_PHYLO_NODE_HPP
#define STS_PARTICLE_DETAIL_PHYLO_NODE_HPP

#include <memory>

#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/TreeTemplateTools.h>

#include "sts/particle/detail/phylo_node_fwd.hpp"
#include "sts/particle/detail/edge_fwd.hpp"
#include "sts/likelihood/detail/online_calculator_fwd.hpp"

namespace sts
{
namespace particle
{

// Implementation
phylo_node::phylo_node(std::shared_ptr<likelihood::online_calculator> calc) : calc(calc) {};
phylo_node::phylo_node(const phylo_node & other) : calc(other.calc)
{
    if(!other.is_leaf()) {
        child1 = std::make_shared<edge>(other.child1->node, other.child1->length);
        child2 = std::make_shared<edge>(other.child2->node, other.child2->length);
    }
}

phylo_node::~phylo_node()
{
    auto p = calc.lock();
    if(p)
        p->unregister_node(this);
}

phylo_node & phylo_node::operator=(const phylo_node & other)
{
    calc = other.calc;
    if(!other.is_leaf()) {
        child1 = std::make_shared<edge>(other.child1->node, other.child1->length);
        child2 = std::make_shared<edge>(other.child1->node, other.child1->length);
    }
    return *this;
}


bool phylo_node::is_leaf() const
{
    return this->child1 == NULL && this->child2 == NULL;
}

void phylo_node::calc_height()
{
    if(is_leaf())
        this->height = 0.0;
    else
        this->height = std::max(child1->node->height + 2 * child1->length, child2->node->height + 2 * child2->length);
}

node phylo_node::of_tree(std::shared_ptr<likelihood::online_calculator> calc, bpp::TreeTemplate<bpp::Node> &tree, int node_number, std::unordered_map<node, std::string>& names)
{
    node n = std::make_shared<phylo_node>(calc);
    if(tree.isLeaf(node_number)) {
        names[n] = tree.getNodeName(node_number);
        return n;
    }
    std::vector<int> children = tree.getSonsId(node_number);
    assert(children.size() == 2);
    n->child1 = edge::of_tree(calc, tree, children[0], names);
    n->child2 = edge::of_tree(calc, tree, children[1], names);
    return n;
}

node phylo_node::of_tree(std::shared_ptr<likelihood::online_calculator> calc, bpp::TreeTemplate<bpp::Node> &tree, std::unordered_map<node, std::string>& names)
{
    return phylo_node::of_tree(calc, tree, tree.getRootId(), names);
}

} // namespace particle
} // namespace sts

#endif // STS_PARTICLE_DETAIL_PHYLO_NODE_HPP

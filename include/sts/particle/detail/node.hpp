/// \file detail/node.hpp

/// \brief phylo_node definition.

#ifndef STS_PARTICLE_DETAIL_PHYLO_NODE_HPP
#define STS_PARTICLE_DETAIL_PHYLO_NODE_HPP

#include <memory>

#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/TreeTemplateTools.h>

#include "sts/particle/detail/node_fwd.hpp"
#include "sts/particle/detail/edge_fwd.hpp"
#include "sts/likelihood/detail/online_calculator_fwd.hpp"

namespace sts
{
namespace particle
{

// Implementation
Node::Node(std::shared_ptr<likelihood::Online_calculator> calc) : calc(calc) {};
Node::Node(const Node & other) : calc(other.calc)
{
    if(!other.is_leaf()) {
        child1 = std::make_shared<Edge>(other.child1->node, other.child1->length);
        child2 = std::make_shared<Edge>(other.child2->node, other.child2->length);
    }
}

Node::~Node()
{
    auto p = calc.lock();
    if(p)
        p->unregister_node(this);
}

Node & Node::operator=(const Node & other)
{
    calc = other.calc;
    if(!other.is_leaf()) {
        child1 = std::make_shared<Edge>(other.child1->node, other.child1->length);
        child2 = std::make_shared<Edge>(other.child1->node, other.child1->length);
    }
    return *this;
}


bool Node::is_leaf() const
{
    return this->child1 == NULL && this->child2 == NULL;
}

node_ptr Node::of_tree(std::shared_ptr<likelihood::Online_calculator> calc, bpp::TreeTemplate<bpp::Node> &tree, int node_number, std::unordered_map<node_ptr, std::string>& names)
{
    node_ptr n = std::make_shared<Node>(calc);
    if(tree.isLeaf(node_number)) {
        names[n] = tree.getNodeName(node_number);
        return n;
    }
    std::vector<int> children = tree.getSonsId(node_number);
    assert(children.size() == 2);
    n->child1 = Edge::of_tree(calc, tree, children[0], names);
    n->child2 = Edge::of_tree(calc, tree, children[1], names);
    return n;
}

node_ptr Node::of_tree(std::shared_ptr<likelihood::Online_calculator> calc, bpp::TreeTemplate<bpp::Node> &tree, std::unordered_map<node_ptr, std::string>& names)
{
    return Node::of_tree(calc, tree, tree.getRootId(), names);
}

} // namespace particle
} // namespace sts

#endif // STS_PARTICLE_DETAIL_PHYLO_NODE_HPP

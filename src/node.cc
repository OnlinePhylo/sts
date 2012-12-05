/// \file node.cc
/// \brief Node implementation

#include "node.h"
#include "edge.h"
#include "online_calculator.h"

namespace sts
{
namespace particle
{

// Implementation
Node::Node(std::shared_ptr<likelihood::Online_calculator> calc) : calc(calc) {}
Node::Node(const Node & other) : calc(other.calc)
{
    if(!other.is_leaf()) {
        child1 = std::make_shared<Edge>(*other.child1);
        child2 = std::make_shared<Edge>(*other.child2);
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
        child1 = std::make_shared<Edge>(*other.child1);
        child2 = std::make_shared<Edge>(*other.child2);
    }
    return *this;
}

/// Is this a leaf node?
bool Node::is_leaf() const
{
    return this->child1 == nullptr && this->child2 == nullptr;
}

Node_ptr Node::of_tree(std::shared_ptr<likelihood::Online_calculator> calc, bpp::TreeTemplate<bpp::Node> &tree, int node_number, std::unordered_map<Node_ptr, std::string>& names)
{
    Node_ptr n = std::make_shared<Node>(calc);
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

Node_ptr Node::of_tree(std::shared_ptr<likelihood::Online_calculator> calc, bpp::TreeTemplate<bpp::Node> &tree, std::unordered_map<Node_ptr, std::string>& names)
{
    return Node::of_tree(calc, tree, tree.getRootId(), names);
}


/// Convenience method
/// \returns the sum of the child edge prior log likelihoods, 0 for leaves.
double Node::edge_prior_log_likelihood() const
{
    if(is_leaf())
        return 0.0;
    return this->child1->prior_log_likelihood + this->child2->prior_log_likelihood;
}

} // namespace particle
} // namespace sts

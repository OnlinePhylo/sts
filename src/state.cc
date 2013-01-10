#include "state.h"
#include "edge.h"
#include "node.h"

namespace sts
{
namespace particle
{

Particle State::of_tree(bpp::TreeTemplate<bpp::Node>& tree,
                        std::unordered_map<sts::particle::Node_ptr, std::string>& names)
{
    Particle p = std::make_shared<State>();
    p->node = Node::of_tree(tree, names);
    if(p->node->is_leaf())
        return p;

    Particle prev = p;
    std::stack<sts::particle::Node_ptr> node_stack;
    node_stack.push(p->node->child1->node);
    node_stack.push(p->node->child2->node);
    while(!node_stack.empty()) {
        Particle cur = std::make_shared<State>();
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

Particle State::of_newick_string(std::string &tree_string,
                                 std::unordered_map<sts::particle::Node_ptr, std::string>& names)
{
    std::unique_ptr<bpp::TreeTemplate<bpp::Node>> tree(bpp::TreeTemplateTools::parenthesisToTree(tree_string));
    Particle node = State::of_tree(*tree, names);
    return node;
}

} // namespace particle
} // namespace sts

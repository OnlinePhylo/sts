/// \file phylo_node_fwd.hpp
/// \brief Forward declarations for a Node.

#ifndef STS_PARTICLE_DETAIL_PHYLO_NODE_FWD_HPP
#define STS_PARTICLE_DETAIL_PHYLO_NODE_FWD_HPP

#include <memory>
#include <string>
#include <unordered_map>
#include <Bpp/Phyl/Node.h>
#include <Bpp/Phyl/TreeTemplate.h>

namespace sts
{

// Circular dependencies
namespace likelihood
{
class online_calculator;
}

namespace particle
{

class edge;

/// \class Node

/// Represents the merge of two trees in a forest.
class Node
{
public:
    explicit Node(std::shared_ptr<likelihood::online_calculator> calc);
    Node(const Node & other);
    ~Node();

    Node & operator=(const Node & other);

    std::shared_ptr<edge> child1;
    std::shared_ptr<edge> child2;

    bool is_leaf() const;

    /// Make a Node from a bpp Tree
    static std::shared_ptr<Node>
    of_tree(std::shared_ptr<likelihood::online_calculator>, bpp::TreeTemplate<bpp::Node> &, std::unordered_map<std::shared_ptr<Node>, std::string>&);

    /// Make a Node from a bpp Tree and node number
    static std::shared_ptr<Node>
    of_tree(std::shared_ptr<likelihood::online_calculator>, bpp::TreeTemplate<bpp::Node> &, int,
            std::unordered_map<std::shared_ptr<Node>, std::string>&);
private:
    std::weak_ptr<likelihood::online_calculator> calc;
};

/// A node in a phylogenetic tree
typedef std::shared_ptr<Node> node_ptr;

}
}

#endif // STS_PARTICLE_DETAIL_PHYLO_NODE_FWD_HPP

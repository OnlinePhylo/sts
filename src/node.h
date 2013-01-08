#ifndef STS_PARTICLE_NODE_H
#define STS_PARTICLE_NODE_H

#include <memory>
#include <string>
#include <unordered_map>
#include <Bpp/Phyl/Node.h>
#include <Bpp/Phyl/TreeTemplate.h>

namespace sts
{

namespace particle
{

class Edge;

/// \class Node
/// \brief Represents the merge of two trees in a forest.
class Node
{
public:
    Node();
    Node(const Node & other);

    Node & operator=(const Node & other);

    std::shared_ptr<Edge> child1;
    std::shared_ptr<Edge> child2;

    bool is_leaf() const;

    /// Make a Node from a bpp Tree
    static std::shared_ptr<Node>
    of_tree(bpp::TreeTemplate<bpp::Node>&,
            std::unordered_map<std::shared_ptr<Node>, std::string>&);

    /// Make a Node from a bpp tree and node number
    static std::shared_ptr<Node>
    of_tree(bpp::TreeTemplate<bpp::Node> &, int,
            std::unordered_map<std::shared_ptr<Node>, std::string>&);

    double edge_prior_log_likelihood() const;
};

/// A node in a phylogenetic tree
typedef std::shared_ptr<Node> Node_ptr;

}
}

#endif // STS_PARTICLE_NODE_H

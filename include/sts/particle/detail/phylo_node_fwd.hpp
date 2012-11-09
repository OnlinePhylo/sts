/// \file phylo_node_fwd.hpp
/// \brief Forward declarations for a phylo_node.

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

/// \class phylo_node

/// Represents the merge of two trees in a forest.
class phylo_node
{
public:
    explicit phylo_node(std::shared_ptr<likelihood::online_calculator> calc);
    phylo_node(const phylo_node & other);
    ~phylo_node();

    phylo_node & operator=(const phylo_node & other);

    std::shared_ptr<edge> child1;
    std::shared_ptr<edge> child2;

    bool is_leaf() const;

    /// Make a phylo_node from a bpp Tree
    static std::shared_ptr<phylo_node>
    of_tree(std::shared_ptr<likelihood::online_calculator>, bpp::TreeTemplate<bpp::Node> &, std::unordered_map<std::shared_ptr<phylo_node>, std::string>&);

    /// Make a phylo_node from a bpp Tree and node number
    static std::shared_ptr<phylo_node>
    of_tree(std::shared_ptr<likelihood::online_calculator>, bpp::TreeTemplate<bpp::Node> &, int,
            std::unordered_map<std::shared_ptr<phylo_node>, std::string>&);
private:
    std::weak_ptr<likelihood::online_calculator> calc;
};

/// A node in a phylogenetic tree
typedef std::shared_ptr<phylo_node> node_ptr;

}
}

#endif // STS_PARTICLE_DETAIL_PHYLO_NODE_FWD_HPP

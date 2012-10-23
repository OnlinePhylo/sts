/// \file phylo_node_fwd.hpp
/// \brief Forward declarations for a phylo_node.

#ifndef STS_PARTICLE_DETAIL_PHYLO_NODE_FWD_HPP
#define STS_PARTICLE_DETAIL_PHYLO_NODE_FWD_HPP

#include <memory>
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
    explicit phylo_node(const phylo_node &other);
    ~phylo_node();

    std::shared_ptr<edge> child1;
    std::shared_ptr<edge> child2;

    // convenience for proposals, height must always increase.
    // In the non-clock case, height is the diameter (2 * distance to closest leaf)
    double height;
    bool is_leaf();

    /// Calculate the height once children have been set
    void calc_height();

    /// Make a phylo_node from a bpp Tree
    static std::shared_ptr<phylo_node>
    of_tree(std::shared_ptr<likelihood::online_calculator>, bpp::TreeTemplate<bpp::Node> &);

    /// Make a phylo_node from a bpp Tree and node number
    static std::shared_ptr<phylo_node>
    of_tree(std::shared_ptr<likelihood::online_calculator>, bpp::TreeTemplate<bpp::Node> &, int);

private:
    std::weak_ptr<likelihood::online_calculator> calc;
};

/// A node in a phylogenetic tree
typedef std::shared_ptr<phylo_node> node;

}
}

#endif // STS_PARTICLE_DETAIL_PHYLO_NODE_FWD_HPP

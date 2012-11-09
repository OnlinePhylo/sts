/// \file edge.cc
/// \brief Edge class
#include "edge.h"
#include "node.h"

namespace sts
{
namespace particle
{
/// Initialize with a node and distance.
///
/// \param node Node on distal size of edge.
/// \param length Branch length.
Edge::Edge(std::shared_ptr<Node> node, double length) : length(length), node(node) {}

/// Initialize with a node only.
///
/// Branch length is initialized at 0.
/// \param node Node on distal size of edge.
Edge::Edge(std::shared_ptr<Node> node) : length(0.0), node(node) {}


std::shared_ptr<Edge> Edge::of_tree(std::shared_ptr<sts::likelihood::Online_calculator> calc, bpp::TreeTemplate<bpp::Node> &tree, int node_number, std::unordered_map<std::shared_ptr<Node>, std::string>& names)
{
    return std::make_shared<Edge>(
               Node::of_tree(calc, tree, node_number, names),
               tree.getDistanceToFather(node_number));
}

} // namespace particle
} // namespace sts

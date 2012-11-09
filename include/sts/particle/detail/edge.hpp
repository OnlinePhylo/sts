/// \file sts/particle/detail/edge.hpp
/// \brief Edge class

#ifndef STS_PARTICLE_DETAIL_EDGE_HPP
#define STS_PARTICLE_DETAIL_EDGE_HPP

#include <memory>

#include "sts/particle/detail/edge_fwd.hpp"
#include "sts/particle/detail/phylo_node_fwd.hpp"
#include "sts/likelihood/detail/online_calculator_fwd.hpp"

namespace sts
{
namespace particle
{
/// Initialize with a node and distance.
///
/// \param node Node on distal size of edge.
/// \param length Branch length.
edge::edge(std::shared_ptr<Node> node, double length) : length(length), node(node) {}

/// Initialize with a node only.
///
/// Branch length is initialized at 0.
/// \param node Node on distal size of edge.
edge::edge(std::shared_ptr<Node> node) : length(0.0), node(node) {}


std::shared_ptr<edge> edge::of_tree(std::shared_ptr<sts::likelihood::online_calculator> calc, bpp::TreeTemplate<bpp::Node> &tree, int node_number, std::unordered_map<std::shared_ptr<Node>, std::string>& names)
{
    return std::make_shared<edge>(
               Node::of_tree(calc, tree, node_number, names),
               tree.getDistanceToFather(node_number));
}

} // namespace particle
} // namespace sts

#endif // STS_PARTICLE_DETAIL_EDGE_HPP

/// \file edge.hpp
/// \brief Edge class

#ifndef STS_PARTICLE_EDGE_HPP
#define STS_PARTICLE_EDGE_HPP

#include <memory>

#include "sts/likelihood/online_calculator_fwd.hpp"
#include "sts/particle/edge_fwd.hpp"

namespace sts
{
namespace particle
{
edge::edge(std::shared_ptr<phylo_node> node, double length) : length(length), node(node) {}

std::shared_ptr<edge> edge::of_tree(std::shared_ptr<sts::likelihood::online_calculator> calc, bpp::TreeTemplate<bpp::Node> &tree, int node_number)
{
    return std::make_shared<edge>(
        phylo_node::of_tree(calc, tree, node_number),
        tree.getDistanceToFather(node_number));
}

} // namespace particle
} // namespace sts

#endif // STS_PARTICLE_EDGE_HPP

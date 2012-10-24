/// \file edge.hpp
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
edge::edge(std::shared_ptr<phylo_node> node, double length) : length(length), node(node) {}

std::shared_ptr<edge> edge::of_tree(std::shared_ptr<sts::likelihood::online_calculator> calc, bpp::TreeTemplate<bpp::Node> &tree, int node_number, std::unordered_map<std::shared_ptr<phylo_node>, std::string>& names)
{
    return std::make_shared<edge>(
               phylo_node::of_tree(calc, tree, node_number, names),
               tree.getDistanceToFather(node_number));
}

} // namespace particle
} // namespace sts

#endif // STS_PARTICLE_DETAIL_EDGE_HPP

#ifndef STS_PARTICLE_DETAIL_EDGE_FWD_HPP
#define STS_PARTICLE_DETAIL_EDGE_FWD_HPP

#include <memory>
#include <string>
#include <unordered_map>
#include <Bpp/Phyl/Tree.h>

namespace sts
{
namespace particle
{
class phylo_node;

/// \class edge
/// An edge
class edge
{
public:
    /// Initialize with a node and distance.
    edge(std::shared_ptr<phylo_node>, double);

    double length;
    std::shared_ptr<phylo_node> node;

    /// Make an edge from a bpp Tree and node number
    static std::shared_ptr<edge> of_tree(std::shared_ptr<sts::likelihood::online_calculator>, bpp::TreeTemplate<bpp::Node> &, int, std::unordered_map<std::shared_ptr<phylo_node>, std::string>&);
};
}
}
#endif

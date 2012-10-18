#ifndef STS_PARTICLE_EDGE_FWD_HPP
#define STS_PARTICLE_EDGE_FWD_HPP

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
    static std::shared_ptr< edge > of_tree(std::shared_ptr<sts::likelihood::online_calculator >, bpp::TreeTemplate<bpp::Node> &, int);
};
}
}
#endif

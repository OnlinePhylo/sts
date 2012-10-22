#ifndef STS_PARTICLE_DETAIL_EDGE_FWD_HPP
#define STS_PARTICLE_DETAIL_EDGE_FWD_HPP

#include <memory>
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
    edge(std::shared_ptr<phylo_node>, double);

    edge(std::shared_ptr<phylo_node>);

    double length;
    std::shared_ptr<phylo_node> node;

    /// Make an edge from a bpp Tree and node number
    static std::shared_ptr<edge> of_tree(std::shared_ptr<sts::likelihood::online_calculator>, bpp::TreeTemplate<bpp::Node> &, int);
};
}
}
#endif

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
class Node;

/// \class edge
/// An edge
class edge
{
public:
    edge(std::shared_ptr<Node>, double);

    edge(std::shared_ptr<Node>);

    /// Length of the edge
    double length;

    /// The node below this edge
    std::shared_ptr<Node> node;

    /// Make an edge from a bpp Tree and node number
    static std::shared_ptr<edge> of_tree(std::shared_ptr<sts::likelihood::Online_calculator>, bpp::TreeTemplate<bpp::Node> &, int, std::unordered_map<std::shared_ptr<Node>, std::string>&);
};
}
}
#endif

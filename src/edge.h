#ifndef STS_PARTICLE_DETAIL_EDGE_H
#define STS_PARTICLE_DETAIL_EDGE_H

#include <memory>
#include <string>
#include <unordered_map>
#include <Bpp/Phyl/TreeTemplate.h>

namespace sts
{
namespace likelihood
{
class Online_calculator;
}
namespace particle
{
class Node;

/// \class Edge
/// An edge, which consists of a double length and a Node.
class Edge
{
public:
    Edge(std::shared_ptr<Node>, double, double);
    Edge(std::shared_ptr<Node>, double);
    explicit Edge(std::shared_ptr<Node>);

    /// Length of the edge
    double length;

    /// Prior likelihood of this edge length
    double prior_log_likelihood;

    /// The node below this edge
    std::shared_ptr<Node> node;

    /// Make an edge from a bpp Tree and node number
    static std::shared_ptr<Edge> of_tree(std::shared_ptr<sts::likelihood::Online_calculator>,
                                         bpp::TreeTemplate<bpp::Node> &,
                                         int,
                                         std::unordered_map<std::shared_ptr<Node>, std::string>&);
};
}
}
#endif

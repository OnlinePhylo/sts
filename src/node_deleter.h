/// \file node_deleter.h
#ifndef STS_LIKELIHOOD_NODE_DELETER_H
#define STS_LIKELIHOOD_NODE_DELETER_H

#include <memory>

namespace sts
{
// Forwards
namespace particle { class Node; }
namespace likelihood
{
class Online_calculator;

/// \brief Deleter for sts::particle::Node which unregisters nodes from an Online_calculator
/// during deletion
/// 
/// Sample usage:
/// \code
/// std::shared_ptr<sts::likelihood::Online_calculator> c;
/// sts::particle::Node_ptr p = sts::particle::Node_ptr(new sts::particle::Node(), Node_deleter(c));
/// \endcode
///
/// \related Node
/// \related Node_ptr
struct Node_deleter
{
    Node_deleter() {};
    Node_deleter(const std::shared_ptr<Online_calculator>& c) : calc(c) {};

    /// The online calculator to unregister the node with on deletion
    std::weak_ptr<Online_calculator> calc;
    std::default_delete<sts::particle::Node> d;

    void operator()(sts::particle::Node* node);
};

}
}

#endif // STS_LIKELIHOOD_NODE_DELETER_H

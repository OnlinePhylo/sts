/// \file node_ptr.h
/// \brief Typedef for shared_ptr<Node>

#ifndef STS_NODE_PTR_H
#define STS_NODE_PTR_H

#include <memory>

namespace sts
{
namespace particle
{
class Node;

/// \brief A shared_ptr to a Node.

/// NB: if using Node_ptr in conjunction with an Online_calculator instance, use
/// a Node_deleter, e.g.:
/// \code
/// std::shared_ptr<Online_calculator> calc = ...;
/// Node_ptr n = Node_ptr(new Node(), Node_deleter(calc));
/// \endcode
/// \related sts::particle::Node
/// \related sts::likelihood::Node_deleter
typedef std::shared_ptr<Node> Node_ptr;
} // sts
} // particle
#endif // STS_NODE_PTR_H

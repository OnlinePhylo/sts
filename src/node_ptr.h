#ifndef STS_NODE_PTR_H
#define STS_NODE_PTR_H

#include <memory>

namespace sts
{
namespace particle
{
class Node;
typedef std::shared_ptr<Node> Node_ptr;
} // sts
} // particle
#endif // STS_NODE_PTR_H

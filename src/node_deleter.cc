#include "node_deleter.h"
#include "node.h"
#include "online_calculator.h"

using sts::particle::Node;

namespace sts
{
namespace likelihood
{

void Node_deleter::operator()(Node* node)
{
    // Unregister from the calculator
    if(auto p = calc.lock())
        p->unregister_node(node);
    // ... before deleting the node
    d(node);
}


} // namespace likelihood
} // namespace sts

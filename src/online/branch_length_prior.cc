#include "branch_length_prior.h"
#include <vector>
#include <gsl/gsl_randist.h>

using bpp::Node;
using bpp::TreeTemplate;
using std::vector;
using std::function;

namespace sts { namespace online {

BranchLengthPrior::BranchLengthPrior(function<double(double)> logPriorDensity) :
    logPriorDensity(logPriorDensity)
{}


double BranchLengthPrior::operator()(const TreeTemplate<Node>& tree) const
{
    double logDensity = 0;
    for(const Node* node : tree.getNodes()) {
        if(!node->hasDistanceToFather())
            continue;
        const double d = node->getDistanceToFather();
        logDensity += logPriorDensity(d);
    }
    return logDensity;
}

}} // Namespaces

#include "branch_length_prior.h"
#include <vector>
#include <gsl/gsl_randist.h>

using bpp::Node;
using bpp::TreeTemplate;
using std::vector;
using std::function;

namespace sts { namespace online {

Branch_length_prior::Branch_length_prior(function<double(double)> log_prior_density) :
    log_prior_density(log_prior_density)
{}


double Branch_length_prior::operator()(const TreeTemplate<Node>& tree) const
{
    double log_density = 0;
    for(const Node* node : tree.getNodes()) {
        if(!node->hasDistanceToFather())
            continue;
        const double d = node->getDistanceToFather();
        log_density += log_prior_density(d);
    }
    return log_density;
}

}} // Namespaces

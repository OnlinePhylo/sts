#include "multiplier_smc_move.h"
#include "composite_tree_likelihood.h"
#include "multiplier_proposal.h"
#include "online_util.h"
#include "tree_particle.h"
#include <gsl/gsl_randist.h>
#include <iostream>
#include <cmath>
#include <utility>
#include <stdexcept>

using namespace bpp;

namespace sts { namespace online {

MultiplierSMCMove::MultiplierSMCMove(CompositeTreeLikelihood& calculator,
                                     const double lambda) :
    calculator(calculator),
    lambda(lambda)
{}

void MultiplierSMCMove::operator()(long, smc::particle<TreeParticle>& particle, smc::rng* rng)
{
    // Choose an edge at random
    TreeParticle* value = particle.GetValuePointer();
    std::vector<bpp::Node*> nodes = onlineAvailableEdges(*value->tree);
    size_t idx = rng->UniformDiscrete(0, nodes.size() - 1);

    bpp::Node* n = nodes[idx];
    const double orig_dist = n->getDistanceToFather();

    calculator.initialize(*value->model, *value->rateDist, *value->tree);

    double orig_ll = calculator();

    const Proposal p = positive_real_multiplier(orig_dist, 1e-6, 100.0, lambda, rng);
    n->setDistanceToFather(p.value);
    double new_ll = calculator();

    particle.AddToLogWeight(new_ll - orig_ll - std::log(p.forwardDensity));
}

}}

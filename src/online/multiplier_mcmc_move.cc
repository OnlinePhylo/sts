#include "multiplier_mcmc_move.h"
#include "composite_tree_likelihood.h"
#include "multiplier_proposal.h"
#include "online_util.h"
#include <cmath>
#include <iostream>
#include <utility>
#include <stdexcept>

using namespace bpp;

namespace sts { namespace online {

MultiplierMCMCMove::MultiplierMCMCMove(CompositeTreeLikelihood& calculator,
                                           const double lambda) :
    OnlineMCMCMove(lambda),
    calculator(calculator)
{}

MultiplierMCMCMove::~MultiplierMCMCMove()
{
    // Debug bits
    if(n_attempted > 0) {
        std::clog << "Multiplier_mcmc_move: " << n_accepted << '/' << n_attempted << ": " << acceptanceProbability() << std::endl;
    }
}

int MultiplierMCMCMove::proposeMove(long, smc::particle<TreeParticle>& particle, smc::rng* rng)
{
    // Choose an edge at random
    TreeParticle* value = particle.GetValuePointer();
    return proposeMove(*value, rng);
}
    
int MultiplierMCMCMove::proposeMove(TreeParticle& particle, smc::rng* rng){
    std::vector<bpp::Node*> nodes = onlineAvailableEdges(*particle.tree);
    size_t idx = rng->UniformDiscrete(0, nodes.size() - 1);
    
    bpp::Node* n = nodes[idx];
    const double orig_dist = n->getDistanceToFather();
    
    calculator.initialize(*particle.model, *particle.rateDist, *particle.tree);
    
    double orig_ll = calculator();
    
    const Proposal p = positive_real_multiplier(orig_dist, 1e-6, 100.0, _lambda, rng);
    n->setDistanceToFather(p.value);
    double new_ll = calculator();
    
    particle.logP = new_ll;
    double mh_ratio = std::exp(new_ll + std::log(p.hastingsRatio) - orig_ll);
    if(mh_ratio >= 1.0 || rng->UniformS() < mh_ratio) {
        return 1;
    } else {
        // Rejected
        particle.logP = orig_ll;
        n->setDistanceToFather(orig_dist);
        return 0;
    }
}
    
}}

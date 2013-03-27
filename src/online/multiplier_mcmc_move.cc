#include "multiplier_mcmc_move.h"
#include "beagle_tree_likelihood.h"
#include "multiplier_proposal.h"
#include "online_util.h"
#include <cmath>
#include <iostream>
#include <utility>
#include <stdexcept>

using namespace bpp;

namespace sts { namespace online {

Multiplier_mcmc_move::Multiplier_mcmc_move(Beagle_tree_likelihood& calculator,
                                           const double lambda) :
    calculator(calculator),
    lambda(lambda)
{}

Multiplier_mcmc_move::~Multiplier_mcmc_move()
{
    // Debug bits
    if(n_attempted > 0) {
        std::clog << "Multiplier_mcmc_move: " << n_accepted << '/' << n_attempted << ": " << acceptance_probability() << std::endl;
    }
}

int Multiplier_mcmc_move::propose_move(long time, smc::particle<Tree_particle>& particle, smc::rng* rng)
{
    // Choose an edge at random
    Tree_particle* value = particle.GetValuePointer();
    std::vector<bpp::Node*> nodes = online_available_edges(*value->tree);
    size_t idx = rng->UniformDiscrete(0, nodes.size() - 1);

    bpp::Node* n = nodes[idx];
    const double orig_dist = n->getDistanceToFather();

    calculator.initialize(*value->model, *value->rate_dist, *value->tree);

    double orig_ll = calculator.calculate_log_likelihood();

    const Proposal p = pos_real_multiplier(orig_dist, 1e-6, 100.0, lambda, rng);
    n->setDistanceToFather(p.value);
    double new_ll = calculator.calculate_log_likelihood();

    double mh_ratio = std::exp(new_ll + std::log(p.hastings_ratio) - orig_ll);
    if(mh_ratio >= 1.0 || rng->UniformS() < mh_ratio) {
        particle.AddToLogWeight(std::log(p.hastings_ratio));
        return 1;
    } else {
        // Rejected
        n->setDistanceToFather(orig_dist);
        return 0;
    }
}

}}

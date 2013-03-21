#include "multiplier_mcmc_move.h"
#include "beagle_tree_likelihood.h"
#include "multiplier_proposal.h"
#include <cmath>
#include <iostream>
#include <utility>
#include <stdexcept>

using namespace bpp;

namespace sts { namespace online {

Multiplier_mcmc_move::Multiplier_mcmc_move(Beagle_tree_likelihood& calculator,
                                           const double lambda) :
    calculator(calculator),
    lambda(lambda),
    proposed(0),
    accepted(0)
{}

Multiplier_mcmc_move::~Multiplier_mcmc_move()
{
    // Debug bits
    if(proposed > 0) {
        double rate = ((double)accepted) / (double)(proposed + accepted);
        std::clog << "Multiplier_mcmc_move: " << accepted << '/' << accepted + proposed << ": " << rate << std::endl;
    }
}

int Multiplier_mcmc_move::operator()(long time, smc::particle<Tree_particle>& particle, smc::rng* rng)
{
    ++proposed;
    // Choose an edge at random
    TreeTemplate<bpp::Node>* tree = particle.GetValuePointer()->tree.get();
    std::vector<bpp::Node*> nodes = online_available_edges(*tree);
    size_t idx = rng->UniformDiscrete(0, nodes.size() - 1);

    bpp::Node* n = nodes[idx];
    const double orig_dist = n->getDistanceToFather();

    calculator.initialize(*particle.GetValuePointer()->model,
                          *particle.GetValuePointer()->rate_dist,
                          *particle.GetValuePointer()->tree);
    double orig_ll = calculator.calculate_log_likelihood();

    const Proposal p = pos_real_multiplier(orig_dist, 1e-6, 100.0, lambda, rng);
    n->setDistanceToFather(p.value);
    double new_ll = calculator.calculate_log_likelihood();

    double mh_ratio = std::exp(new_ll + std::log(p.hastings_ratio) - orig_ll);
    if(mh_ratio >= 1.0 || rng->UniformS() < mh_ratio) {
        ++accepted;
        particle.AddToLogWeight(std::log(p.hastings_ratio));
        return 1;
    } else {
        // Rejected
        n->setDistanceToFather(orig_dist);
        return 0;
    }
}

}}

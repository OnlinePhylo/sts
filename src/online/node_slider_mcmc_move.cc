#include "node_slider_mcmc_move.h"
#include "beagle_tree_likelihood.h"
#include "multiplier_proposal.h"
#include "online_util.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <utility>
#include <unordered_set>

using namespace bpp;

namespace sts { namespace online {

Node_slider_mcmc_move::Node_slider_mcmc_move(Beagle_tree_likelihood& calculator,
                                             const double lambda) :
    calculator(calculator),
    lambda(lambda)
{}

Node_slider_mcmc_move::~Node_slider_mcmc_move()
{
    // Debug bits
    if(n_attempted > 0) {
        std::clog << "Node_slider_mcmc_move: " << n_accepted << '/' << n_attempted << ": " << acceptance_probability() << std::endl;
    }
}

int Node_slider_mcmc_move::propose_move(long time, smc::particle<Tree_particle>& particle, smc::rng* rng)
{
    // Choose an edge at random
    TreeTemplate<bpp::Node>* tree = particle.GetValuePointer()->tree.get();
    std::vector<bpp::Node*> nodes = online_available_edges(*tree);

    // restrict to nodes with a parent in the available set
    const std::unordered_set<Node*> node_set(nodes.begin(), nodes.end());
    auto father_unavailable = [&node_set](Node* node) {
        return node_set.find(node->getFather()) == node_set.end();
    };
    nodes.erase(std::remove_if(nodes.begin(), nodes.end(), father_unavailable), nodes.end());

    size_t idx = rng->UniformDiscrete(0, nodes.size() - 1);
    Node* n = nodes[idx];

    Node* father = n->getFather();

    const double orig_dist = n->getDistanceToFather() + father->getDistanceToFather();

    calculator.initialize(*particle.GetValuePointer()->model,
                          *particle.GetValuePointer()->rate_dist,
                          *tree);
    double orig_ll = calculator.calculate_log_likelihood();

    const Proposal p = pos_real_multiplier(orig_dist, 1e-6, 100.0, lambda, rng);
    const double d = rng->UniformS() * p.value;

    n->setDistanceToFather(d);
    father->setDistanceToFather(p.value - d);

    double new_ll = calculator.calculate_log_likelihood();

    double mh_ratio = std::exp(new_ll + std::log(p.hastings_ratio) - orig_ll);
    if(mh_ratio >= 1.0 || rng->UniformS() < mh_ratio) {
        particle.AddToLogWeight(-std::log(p.hastings_ratio));
        return 1;
    } else {
        // Rejected
        n->setDistanceToFather(orig_dist);
        return 0;
    }
}

}}

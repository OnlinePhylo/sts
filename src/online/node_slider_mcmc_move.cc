#include "node_slider_mcmc_move.h"
#include "composite_tree_likelihood.h"
#include "multiplier_proposal.h"
#include "online_util.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <utility>
#include <unordered_set>

using namespace bpp;

namespace sts { namespace online {

NodeSliderMCMCMove::NodeSliderMCMCMove(CompositeTreeLikelihood& calculator,
                                       const double lambda) :
    OnlineMCMCMove(lambda),
    calculator(calculator)
{}

NodeSliderMCMCMove::~NodeSliderMCMCMove()
{
    // Debug bits
    if(n_attempted > 0) {
        std::clog << "Node_slider_mcmc_move: " << n_accepted << '/' << n_attempted << ": " << acceptanceProbability() << std::endl;
    }
}

int NodeSliderMCMCMove::proposeMove(long, smc::particle<TreeParticle>& particle, smc::rng* rng)
{
    // Choose an edge at random
    TreeParticle* value = particle.GetValuePointer();
    TreeTemplate<bpp::Node>* tree = particle.GetValuePointer()->tree.get();
    std::vector<bpp::Node*> nodes = onlineAvailableEdges(*tree);

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

    const double orig_n_dist = n->getDistanceToFather();
    const double orig_father_dist = father->getDistanceToFather();

    calculator.initialize(*particle.GetValuePointer()->model,
                          *particle.GetValuePointer()->rateDist,
                          *tree);
    double orig_ll = calculator();

    const Proposal p = positive_real_multiplier(orig_dist, 1e-6, 100.0, _lambda, rng);
    const double d = rng->UniformS() * p.value;

    n->setDistanceToFather(d);
    father->setDistanceToFather(p.value - d);

    double new_ll = calculator();
    value->logP = new_ll;

    double mh_ratio = std::exp(new_ll + std::log(p.hastingsRatio) - orig_ll);
    if(mh_ratio >= 1.0 || rng->UniformS() < mh_ratio) {
        return 1;
    } else {
        // Rejected
        value->logP = orig_ll;
        n->setDistanceToFather(orig_n_dist);
        father->setDistanceToFather(orig_father_dist);
        return 0;
    }
}

}}

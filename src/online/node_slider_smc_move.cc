#include "node_slider_smc_move.h"
#include "composite_tree_likelihood.h"
#include "multiplier_proposal.h"
#include "tree_particle.h"
#include "online_util.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <utility>
#include <unordered_set>

using namespace bpp;

namespace sts { namespace online {

NodeSliderSMCMove::NodeSliderSMCMove(CompositeTreeLikelihood& calculator,
                                     const double lambda) :
    calculator(calculator),
    lambda(lambda)
{}

void NodeSliderSMCMove::operator()(long, smc::particle<TreeParticle>& particle, smc::rng* rng)
{
    // Choose an edge at random
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

    calculator.initialize(*particle.GetValuePointer()->model,
                          *particle.GetValuePointer()->rateDist,
                          *tree);
    double orig_ll = calculator();

    const Proposal p = positive_real_multiplier(orig_dist, 1e-6, 1000.0, lambda, rng);
    const double d = rng->UniformS() * p.value;

    n->setDistanceToFather(d);
    father->setDistanceToFather(p.value - d);

    double new_ll = calculator();

    particle.AddToLogWeight(new_ll - orig_ll - p.forwardDensity);
}

}}

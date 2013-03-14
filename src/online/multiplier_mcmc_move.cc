#include "multiplier_mcmc_move.h"
#include "beagle_tree_likelihood.h"
#include <cmath>
#include <utility>

using namespace bpp;

namespace sts { namespace online {

Multiplier_mcmc_move::Multiplier_mcmc_move(Beagle_tree_likelihood& calculator,
                                           const double lambda) :
    calculator(calculator),
    lambda(lambda)
{}

const double DEFAULT_TUNE = 2.0 * std::log(1.2);

typedef std::pair<double, double> proposal_hastings;

struct Proposal
{
    double value;
    double hastings_ratio;
};

Proposal pos_real_multiplier(const double value, const double min_value, const double max_value, const double tuning, smc::rng* rng)
{
    const double random = rng->UniformS();
    const double factor = std::exp(tuning * (random - 0.5));
    double new_value = factor * value;
    bool valid = false;
    do {
        if(new_value < min_value)
            new_value = min_value * min_value / new_value;
        else if(new_value > max_value)
            new_value = max_value * max_value / new_value;
        else
            valid = true;
    } while (!valid);

    return {new_value, new_value / value};
}

int Multiplier_mcmc_move::operator()(long time, smc::particle<Tree_particle>& particle, smc::rng* rng)
{
    ++proposed;
    // Choose an edge at random
    TreeTemplate<Node>* tree = particle.GetValuePointer()->tree.get();
    std::vector<Node*> nodes = tree->getNodes();
    size_t idx = rng->UniformDiscrete(0, nodes.size() - 2);
    if(nodes[idx] == tree->getRootNode())
        idx++;

    Node* n = nodes[idx];
    const double orig_dist = n->getDistanceToFather();

    calculator.load_rate_distribution(*particle.GetValuePointer()->rate_dist);
    calculator.load_substitution_model(*particle.GetValuePointer()->model);
    double orig_ll = calculator.calculate_log_likelihood(*tree);
    const Proposal p = pos_real_multiplier(orig_dist, 1e-6, 100.0, lambda, rng);
    n->setDistanceToFather(p.value);
    double new_ll = calculator.calculate_log_likelihood(*tree);

    double met = new_ll + std::log(p.hastings_ratio) - orig_ll;
    if(met >= 1.0 || rng->UniformS() < met) {
        ++accepted;
        return 1;
    } else {
        // Rejected
        n->setDistanceToFather(orig_dist);
        return 0;
    }
}

}}

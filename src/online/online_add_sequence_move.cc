#include "online_add_sequence_move.h"
#include "tree_particle.h"
#include "beagle_tree_likelihood.h"
#include "likelihood_vector.h"
#include "gsl.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <memory>

#include <gsl/gsl_randist.h>

using namespace std;
using namespace bpp;

namespace sts { namespace online {

double minimize(std::function<double(double)> fn,
                double raw_start,
                double left,
                double right,
                const size_t max_iters=5)
{
    size_t iter = 0;

    double lefty = fn(left);
    double righty = fn(right);
    double start = raw_start;
    double val;
    double min_x = lefty < righty ? left : right;
    double min_y = std::min(righty, lefty);

    for(iter = 0; iter < max_iters; iter++) {
        val = fn(start);
        if(fn(start) < min_y)
            return sts::gsl::minimize(fn, start, left, right, max_iters - iter);

        if(std::abs(start - min_x) < 1e-3)
            return start;
        start = (start + min_x) / 2;
    }

    return start;

}

Online_add_sequence_move::Online_add_sequence_move(Beagle_tree_likelihood& calculator,
                                                   const vector<string>& taxa_to_add) :
    calculator(calculator),
    taxa_to_add(taxa_to_add)
{ }


int Online_add_sequence_move::operator()(long time, smc::particle<Tree_particle>& particle, smc::rng* rng) const
{
    Tree_particle* value = particle.GetValuePointer();
    unique_ptr<TreeTemplate<bpp::Node>>& tree = value->tree;

    const size_t orig_n_leaves = tree->getNumberOfLeaves(),
                 orig_n_nodes = tree->getNumberOfNodes();


    assert(time - 1 >= 0);
    size_t i = time - 1;
    assert(i < taxa_to_add.size());

    // Replace node `n` in the tree with a new node containing as children `n` and `new_node`
    // Attach a new leaf, in the following configuration
    //
    //              father
    //   /          o
    //   |          | d - dist_bl
    //   |          |
    // d | new_node o-------o new_leaf
    //   |          |
    //   |          | dist_bl
    //   \          o
    //              n

    calculator.initialize(*value->model, *value->rate_dist, *tree);


    // Choose an edge for insertion
    // This is a guided move - we calculate the likelihood at the center of each edge
    // TODO: Do we need special handling for the root?
    Likelihood_vector partials = calculator.get_leaf_partials(taxa_to_add[i]);
    const vector<Beagle_tree_likelihood::Node_partials> np = calculator.get_mid_edge_partials();
    vector<double> products;
    products.reserve(np.size());
    for(const auto& i : np) {
        double product = partials.log_dot(i.second);
        products.push_back(product);
    }

    // Subtract the max LL to avoid underflows
    double max_ll = *std::max_element(products.begin(), products.end());
    for(double& p : products)
        p = std::exp(p - max_ll);

    unique_ptr<unsigned[]> indexes(new unsigned[products.size()]);

    // Select an edge
    rng->Multinomial(1, np.size(), products.data(), indexes.get());

    unsigned idx = 0;
    for(size_t i = 0; i < products.size(); i++) {
        if(indexes[i] == 1) {
            idx = i;
            break;
        }
    }

    // Proposal density
    particle.AddToLogWeight(-std::log(products[idx]));

    // New internal node, new leaf
    Node* new_node = new Node();
    Node* new_leaf = new Node(taxa_to_add[i]);
    new_node->addSon(new_leaf);
    new_leaf->setDistanceToFather(rng->Exponential(1.0));

    Node* n = tree->getNode(np[idx].first->getId());
    assert(n->hasFather());
    Node* father = n->getFather();

    // Uniform distribution on attachment location
    const double d = n->getDistanceToFather();
    const double dist_bl = d * rng->Uniform(0.0, 1.0);

    // Swap `new_node` in for `n`
    // Note: use {add,remove}Son, rather than {remove,set}Father -
    // latter functions do not update parent sons list.
    father->addSon(new_node);
    father->removeSon(n);
    new_node->addSon(n);

    // Attachment branch lengths
    new_node->setDistanceToFather(d - dist_bl);
    n->setDistanceToFather(dist_bl);

    assert(!tree->isMultifurcating());
    assert(tree->isRooted());
    assert(new_node->getNumberOfSons() == 2);
    assert(!new_node->isLeaf());
    assert(new_leaf->getNumberOfSons() == 0);
    assert(new_leaf->isLeaf());
    assert(tree->getNumberOfLeaves() == orig_n_leaves + 1);
    assert(tree->getNumberOfNodes() == orig_n_nodes + 2);

    // TODO: Proposal density
    // TODO: Do we need a backward proposal density? Given an (arbitrary) root on an unrooted tree, the (unrooted) edge
    // containing the root will have 2 potential attachment points.
    // Implemented here.
    if(new_node->getFather() == tree->getRootNode())
        particle.AddToLogWeight(std::log(0.5));

    // Calculate new LL
    calculator.initialize(*value->model, *value->rate_dist, *value->tree);

    // A little fitting
    //auto pend_bl_fn = [&](double pend_bl) -> double {
        //double orig_bl = new_leaf->getDistanceToFather();
        //new_leaf->setDistanceToFather(pend_bl);
        //double ll = -calculator.calculate_log_likelihood(*tree);
        //new_leaf->setDistanceToFather(orig_bl);
        //return ll;
    //};
    //double best_pend = minimize(pend_bl_fn, 1.0, 1e-7, 2.0);

    // Propose a pendant branch length from an exponential distribution around best_pend
    new_leaf->setDistanceToFather(rng->Exponential(0.1));
    //particle.AddToLogWeight(-std::log(gsl_ran_exponential_pdf(new_leaf->getDistanceToFather(), best_pend)));

    auto dist_bl_fn = [&](double dist_bl) -> double {
        n->setDistanceToFather(dist_bl);
        new_node->setDistanceToFather(d - dist_bl);
        double ll = -calculator.calculate_log_likelihood();
        return ll;
    };
    const double best_dist = minimize(dist_bl_fn, d/2., 1e-6, d);

    // Propose a distal branch length
    const double new_distal = std::min(rng->Exponential(best_dist), d);
    new_node->setDistanceToFather(d - new_distal);
    n->setDistanceToFather(new_distal);
    //particle.AddToLogWeight(-std::log(gsl_ran_exponential_pdf(new_distal, best_dist)));

    const double log_like = calculator.calculate_log_likelihood();
    particle.AddToLogWeight(log_like);

    return 0;
}

}} // namespaces

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


/// Choose an edge for insertion
/// This is a guided move - we calculate the likelihood with the sequence inserted at the middle of each
/// edge, then select an edge by sampling from the multinomial distribution weighted by the edge-likelihoods.
pair<Node*, double> Online_add_sequence_move::choose_edge(TreeTemplate<Node>& tree, const std::string& leaf_name, smc::rng* rng)
{
    // First, calculate the products
    Likelihood_vector partials = calculator.get_leaf_partials(leaf_name);
    const vector<Beagle_tree_likelihood::Node_partials> np = calculator.get_mid_edge_partials();
    vector<double> edge_log_likes;
    edge_log_likes.reserve(np.size());
    for(const auto& i : np) {
        double edge_log_like = partials.log_dot(i.second);
        edge_log_likes.push_back(edge_log_like);
    }

    // Find & subtract the max LL to avoid underflow, exponentiate
    double max_ll = *std::max_element(edge_log_likes.begin(), edge_log_likes.end());
    vector<double> edge_likes(edge_log_likes.size());
    std::transform(edge_log_likes.begin(), edge_log_likes.end(), edge_likes.begin(),
                   [&max_ll](double p) { return std::exp(p - max_ll); });

    // Select an edge
    std::vector<unsigned> indexes(edge_likes.size());
    rng->Multinomial(1, np.size(), edge_likes.data(), indexes.data());
    // Only one value should be selected - find it
    auto positive = [](const unsigned x) { return x > 0; };
    std::vector<unsigned>::const_iterator it = std::find_if(indexes.begin(), indexes.end(),
                                                            positive);
    assert(it != indexes.end());
    const size_t idx = it - indexes.begin();

    Node* n = tree.getNode(np[idx].first->getId());
    return pair<Node*,double>(n, edge_log_likes[idx]);
}



int Online_add_sequence_move::operator()(long time, smc::particle<Tree_particle>& particle, smc::rng* rng)
{
    Tree_particle* value = particle.GetValuePointer();
    unique_ptr<TreeTemplate<bpp::Node>>& tree = value->tree;

    const size_t orig_n_leaves = tree->getNumberOfLeaves(),
                 orig_n_nodes = tree->getNumberOfNodes();


    assert(time - 1 >= 0);
    const size_t i = time - 1;
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

    pair<Node*,double> edge_lnp = choose_edge(*tree, taxa_to_add[i], rng);

    // Subtract proposal density (this is q(s_{r-1} \rightarrow s_r))
    particle.AddToLogWeight(-edge_lnp.second);

    // New internal node, new leaf
    Node* new_node = new Node();
    Node* new_leaf = new Node(taxa_to_add[i]);
    new_node->addSon(new_leaf);
    new_leaf->setDistanceToFather(rng->Exponential(1.0));

    Node* n = edge_lnp.first;
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

    // Verify some postconditions
    assert(!tree->isMultifurcating());
    assert(tree->isRooted());
    assert(new_node->getNumberOfSons() == 2);
    assert(!new_node->isLeaf());
    assert(new_leaf->getNumberOfSons() == 0);
    assert(new_leaf->isLeaf());
    assert(tree->getNumberOfLeaves() == orig_n_leaves + 1);
    assert(tree->getNumberOfNodes() == orig_n_nodes + 2);

    // Calculate new LL - need to re-initialize since nodes have been added
    // TODO: Should nodes be allocated dynamically?
    calculator.initialize(*value->model, *value->rate_dist, *value->tree);

    // Propose a pendant branch length from an exponential distribution around best_pend
    new_leaf->setDistanceToFather(rng->Exponential(0.1));
    //particle.AddToLogWeight(-std::log(gsl_ran_exponential_pdf(new_leaf->getDistanceToFather(), best_pend)));

    const double log_like = calculator.calculate_log_likelihood();
    particle.AddToLogWeight(log_like);

    return 0;
}

}} // namespaces

#include "online_add_sequence_move.h"
#include "tree_particle.h"
#include "beagle_tree_likelihood.h"
#include "composite_tree_likelihood.h"
#include "likelihood_vector.h"
#include "gsl.h"
#include "util.h"

#include <algorithm>
#include <iterator>
#include <cassert>
#include <cmath>
#include <memory>

#include <gsl/gsl_randist.h>

using namespace std;
using namespace bpp;
using sts::util::beagle_check;

namespace sts { namespace online {

const static double TOLERANCE = 1e-3;

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
        if(val < min_y)
            return sts::gsl::minimize(fn, start, left, right, max_iters - iter);

        if(std::abs(start - min_x) < TOLERANCE)
            return start;
        start = (start + min_x) / 2;
    }

    return start;

}

struct TripodOptimizer
{
    int beagleInstance,
        distalBuffer,
        proximalBuffer,
        leafBuffer,
        scratch1,
        scratch2;
    double d;

    /// Optimize distal branch length, keeping pendant fixed
    double optimize_distal(const double distal_start, const double pendant, size_t max_iters=10)
    {
        auto fn = [&](double distal) {
            return - log_like(distal, pendant, true);
        };
        return minimize(fn, distal_start, 0, d, max_iters);
    }

    /// Optimize pendant branch length, keeping distal fixed
    double optimize_pendant(const double distal, const double pendant_start, size_t max_iters=10)
    {
        auto fn = [&](double pendant) {
            return -log_like(distal, pendant, false);
        };

        return minimize(fn, pendant_start, 0, 2.0, max_iters);
    }

    double log_like(const double distal, const double pendant, const bool distal_changed=true)
    {
        std::vector<BeagleOperation> operations;
        std::vector<double> branch_lengths;
        std::vector<int> node_indices;

        // If distal changed, update partial
        if(distal_changed) {
            operations.push_back(BeagleOperation({scratch1,
                                                 BEAGLE_OP_NONE,
                                                 BEAGLE_OP_NONE,
                                                 distalBuffer,
                                                 distalBuffer,
                                                 proximalBuffer,
                                                 proximalBuffer}));
            branch_lengths.push_back(distal);
            node_indices.push_back(distalBuffer);
            branch_lengths.push_back(d - distal);
            node_indices.push_back(proximalBuffer);
        }
        // Always update root partials
        operations.push_back(BeagleOperation({scratch2,
                                              BEAGLE_OP_NONE,
                                              BEAGLE_OP_NONE,
                                              scratch1,
                                              scratch1,
                                              leafBuffer,
                                              leafBuffer}));
        branch_lengths.push_back(0);
        node_indices.push_back(scratch1);
        branch_lengths.push_back(pendant);
        node_indices.push_back(leafBuffer);

        // Usual thing
        beagle_check(beagleUpdateTransitionMatrices(beagleInstance,
                                                    0,
                                                    node_indices.data(),
                                                    NULL,
                                                    NULL,
                                                    branch_lengths.data(),
                                                    node_indices.size()));
        beagle_check(beagleUpdatePartials(beagleInstance, operations.data(), operations.size(), scratch2));

        std::vector<int> scale_indices(operations.size());
        for(size_t i = 0; i < operations.size(); i++)
            scale_indices[i] = operations[i].destinationPartials;

        beagle_check(beagleAccumulateScaleFactors(beagleInstance, scale_indices.data(), scale_indices.size(),
scratch2));
        const int category_weight_index = 0;
        const int state_frequency_index = 0;
        double log_likelihood;
        beagle_check(beagleCalculateRootLogLikelihoods(beagleInstance,
                                                       &scratch2,
                                                       &category_weight_index,
                                                       &state_frequency_index,
                                                       &scratch2,
                                                       1,
                                                       &log_likelihood));
        return log_likelihood;
    }
};

OnlineAddSequenceMove::OnlineAddSequenceMove(CompositeTreeLikelihood& calculator,
                                             const vector<string>& taxaToAdd) :
    calculator(calculator),
    taxaToAdd(std::begin(taxaToAdd), std::end(taxaToAdd)),
    lastTime(-1)
{ }

/// Choose an edge for insertion
/// This is a guided move - we calculate the likelihood with the sequence inserted at the middle of each
/// edge, then select an edge by sampling from the multinomial distribution weighted by the edge-likelihoods.
pair<Node*, double> OnlineAddSequenceMove::chooseEdge(TreeTemplate<Node>& tree, const std::string& leaf_name,
                                                       smc::rng* rng)
{
    // First, calculate the products
    const int leafBuffer = calculator.calculator()->getLeafBuffer(leaf_name);
    const vector<BeagleTreeLikelihood::NodePartials> np = calculator.calculator()->getMidEdgePartials();
    vector<double> edge_log_likes;
    edge_log_likes.reserve(np.size());
    for(const auto& i : np) {
        const double edgeLogLike = calculator.calculator()->logDot(i.second, leafBuffer);
        edge_log_likes.push_back(edgeLogLike);
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
    return pair<Node*,double>(n, edge_log_likes[idx] - max_ll);
}


/// Propose branch-lengths around ML value
AttachmentLocation OnlineAddSequenceMove::proposeBranchLengths(const Node* insertEdge,
                                                               const std::string& newLeafName)
{
    const double d = insertEdge->getDistanceToFather();

    BeagleTreeLikelihood& btl = *calculator.calculator();

    assert(btl.freeBufferCount() >= 2);
    BeagleBuffer b1 = btl.borrowBuffer(), b2 = btl.borrowBuffer();
    // Initialize
    TripodOptimizer optim;
    optim.beagleInstance = calculator.calculator()->beagleInstance();
    optim.scratch1 = b1.value();
    optim.scratch2 = b2.value();
    optim.distalBuffer = calculator.calculator()->getDistalBuffer(insertEdge);
    optim.proximalBuffer = calculator.calculator()->getProximalBuffer(insertEdge);
    optim.leafBuffer = calculator.calculator()->getLeafBuffer(newLeafName);

    optim.d = d;

    double pendant = 1e-8;
    double distal = d / 2;

    // Optimize distal, pendant up to 5 times
    for(size_t i = 0; i < 5; i++) {
        const double new_distal = optim.optimize_distal(distal, pendant);
        if(std::abs(new_distal - distal) < TOLERANCE)
            break;

        const double new_pendant = optim.optimize_pendant(new_distal, pendant);
        if(std::abs(new_pendant - pendant) < TOLERANCE)
            break;

        pendant = new_pendant;
        distal = new_distal;
    }
    return AttachmentLocation{distal, pendant};
}

void OnlineAddSequenceMove::operator()(long time, smc::particle<TreeParticle>& particle, smc::rng* rng)
{
    if(time != lastTime && lastTime >= 0)
        taxaToAdd.pop_front();
    lastTime = time;

    if(taxaToAdd.empty()) {
        assert(0 && "No more sequences to add");
    }

    TreeParticle* value = particle.GetValuePointer();
    unique_ptr<TreeTemplate<bpp::Node>>& tree = value->tree;

    const size_t orig_n_leaves = tree->getNumberOfLeaves(),
                 orig_n_nodes = tree->getNumberOfNodes();

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

    calculator.initialize(*value->model, *value->rateDist, *tree);

    // Calculate root log-likelihood of original tree
    // \gamma*(s_{r-1,k}) from PhyloSMC eqn 2
    const double orig_ll = calculator();
    particle.AddToLogWeight(-orig_ll);

    pair<Node*,double> edge_lnp = chooseEdge(*tree, taxaToAdd.front(), rng);

    // Subtract proposal density (this is q(s_{r-1} \rightarrow s_r))
    particle.AddToLogWeight(-edge_lnp.second);

    // New internal node, new leaf
    Node* new_node = new Node();
    Node* new_leaf = new Node(taxaToAdd.front());
    new_node->addSon(new_leaf);

    Node* n = edge_lnp.first;
    assert(n->hasFather());
    Node* father = n->getFather();

    // branch lengths
    AttachmentLocation ml_bls = proposeBranchLengths(n, new_leaf->getName());

    const double d = n->getDistanceToFather();
    double dist_bl = -1;

    // Handle very small branch lengths - attach with distal BL of 0
    if(d < 1e-8)
        dist_bl = 0;
    else {
        do {
            dist_bl = rng->NormalTruncated(ml_bls.distal_bl, d / 4, 0.0);
        } while(dist_bl < 0 || dist_bl > d);
    }

    assert(!std::isnan(dist_bl));

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
    calculator.initialize(*value->model, *value->rateDist, *value->tree);

    // Propose a pendant branch length from an exponential distribution around best_pend
    new_leaf->setDistanceToFather(rng->Exponential(ml_bls.pendant_bl));
    //particle.AddToLogWeight(-std::log(gsl_ran_exponential_pdf(new_leaf->getDistanceToFather(), best_pend)));

    const double log_like = calculator();
    particle.AddToLogWeight(log_like);
}

}} // namespaces

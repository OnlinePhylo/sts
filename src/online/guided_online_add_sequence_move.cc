#include "guided_online_add_sequence_move.h"
#include "tree_particle.h"
#include "beagle_tree_likelihood.h"
#include "composite_tree_likelihood.h"
#include "likelihood_vector.h"
#include "gsl.h"
#include "util.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iterator>
#include <limits>

#include <gsl/gsl_randist.h>

using namespace std;
using namespace bpp;
using sts::util::beagle_check;

namespace sts { namespace online {

const static double TOLERANCE = 1e-3;

double minimize(std::function<double(double)> fn,
                double rawStart,
                double left,
                double right,
                const size_t maxIters=5)
{
    size_t iter = 0;

    double lefty = fn(left);
    double righty = fn(right);
    double start = rawStart;
    double val;
    double min_x = lefty < righty ? left : right;
    double min_y = std::min(righty, lefty);

    for(iter = 0; iter < maxIters; iter++) {
        val = fn(start);
        if(val < min_y)
            return sts::gsl::minimize(fn, start, left, right, maxIters - iter);

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
    double optimizeDistal(const double distal_start, const double pendant, size_t max_iters=10)
    {
        auto fn = [&](double distal) {
            return - log_like(distal, pendant, true);
        };
        return minimize(fn, distal_start, 0, d, max_iters);
    }

    /// Optimize pendant branch length, keeping distal fixed
    double optimizePendant(const double distal, const double pendant_start, size_t max_iters=10)
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
        const int categoryWeightIdx = 0;
        const int stateFreqIdx = 0;
        double logLike;
        beagle_check(beagleCalculateRootLogLikelihoods(beagleInstance,
                                                       &scratch2,
                                                       &categoryWeightIdx,
                                                       &stateFreqIdx,
                                                       &scratch2,
                                                       1,
                                                       &logLike));
        return logLike;
    }
};

GuidedOnlineAddSequenceMove::GuidedOnlineAddSequenceMove(CompositeTreeLikelihood& calculator,
                                                         const vector<string>& taxaToAdd,
                                                         const vector<double>& proposePendantBranchLengths) :
    OnlineAddSequenceMove(calculator, taxaToAdd),
    proposePendantBranchLengths(proposePendantBranchLengths)
{
    assert(!proposePendantBranchLengths.empty() && "No proposal branch lengths!");
}

/// Choose an edge for insertion
/// This is a guided move - we calculate the likelihood with the sequence inserted at the middle of each
/// edge, then select an edge by sampling from the multinomial distribution weighted by the edge-likelihoods.
///
/// Likelihoods are calculated for each branch length in #proposePendantBranchLengths. The likelihoods from the branch
/// length with the highest / overall likelihood are used for the proposal.
pair<Node*, double> GuidedOnlineAddSequenceMove::chooseEdge(TreeTemplate<Node>& tree, const std::string& leaf_name,
                                                            smc::rng* rng)
{
    // First, calculate the products
    const int leafBuffer = calculator.calculator()->getLeafBuffer(leaf_name);
    const vector<BeagleTreeLikelihood::NodePartials> np = calculator.calculator()->getMidEdgePartials();
    vector<double> edge_log_likes(np.size(), -std::numeric_limits<double>::max());
    for(const double d : proposePendantBranchLengths) {
        for(size_t i = 0; i < np.size(); i++) {
            const double edgeLogLike = calculator.calculator()->logDot(np[i].second, leafBuffer, d);
            edge_log_likes[i] = std::max(edgeLogLike, edge_log_likes[i]);
        }
    }

    const double bestLogLike = *std::max_element(edge_log_likes.begin(), edge_log_likes.end());
    vector<double> edge_likes(edge_log_likes.size());
    std::transform(edge_log_likes.begin(), edge_log_likes.end(), edge_likes.begin(),
                   [&bestLogLike](double p) { return std::exp(p - bestLogLike); });

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
    return pair<Node*,double>(n, edge_log_likes[idx] - bestLogLike);
}

/// Propose branch-lengths around ML value
void GuidedOnlineAddSequenceMove::optimizeBranchLengths(const Node* insertEdge,
                                                        const std::string& newLeafName,
                                                        double& distalBranchLength,
                                                        double& pendantBranchLength)
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
        const double newDistal = optim.optimizeDistal(distal, pendant);
        if(std::abs(newDistal - distal) < TOLERANCE)
            break;

        const double newPendant = optim.optimizePendant(newDistal, pendant);
        if(std::abs(newPendant - pendant) < TOLERANCE)
            break;

        pendant = newPendant;
        distal = newDistal;
    }

    distalBranchLength = distal;
    pendantBranchLength = pendant;
}

AttachmentProposal GuidedOnlineAddSequenceMove::propose(const std::string& leafName, smc::particle<TreeParticle>& particle, smc::rng* rng)
{
    TreeParticle* value = particle.GetValuePointer();
    unique_ptr<TreeTemplate<bpp::Node>>& tree = value->tree;

    // Replace node `n` in the tree with a new node containing as children `n` and `new_node`
    // Attach a new leaf, in the following configuration
    //
    //              father
    //   /          o
    //   |          | d - distal
    //   |          |
    // d | new_node o-------o new_leaf
    //   |          |
    //   |          | distal
    //   \          o
    //              n

    Node* n = nullptr;
    double edgeLogDensity;
    // branch lengths
    std::tie(n, edgeLogDensity) = chooseEdge(*tree, leafName, rng);
    double mlDistal, mlPendant;
    optimizeBranchLengths(n, leafName, mlDistal, mlPendant);

    const double d = n->getDistanceToFather();
    double distal = -1;

    // Handle very small branch lengths - attach with distal BL of 0
    if(d < 1e-8)
        distal = 0;
    else {
        do {
            distal = rng->NormalTruncated(mlDistal, d / 4, 0.0);
        } while(distal < 0 || distal > d);
    }
    assert(!std::isnan(distal));

    const double distalLogDensity = std::log(gsl_ran_gaussian_pdf(distal - mlDistal, d / 4));
    assert(!std::isnan(distalLogDensity));
    const double pendantBranchLength = rng->Exponential(mlPendant);
    const double pendantLogDensity = std::log(gsl_ran_exponential_pdf(pendantBranchLength, mlPendant));
    assert(!std::isnan(pendantLogDensity));
    return AttachmentProposal { n, edgeLogDensity, distal, distalLogDensity, pendantBranchLength, pendantLogDensity, mlDistal, mlPendant };
}

}} // namespaces

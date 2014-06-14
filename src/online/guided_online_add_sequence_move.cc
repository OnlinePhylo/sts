#include "guided_online_add_sequence_move.h"
#include "tree_particle.h"
#include "beagle_tree_likelihood.h"
#include "composite_tree_likelihood.h"
#include "likelihood_vector.h"
#include "tripod_optimizer.h"
#include "util.h"
#include "online_util.h"
#include "log_tricks.h"
#include "weighted_selector.h"

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


/// \brief Discretize the tree into \c n evenly spaced points
std::vector<BeagleTreeLikelihood::AttachmentLocation> discretizeTree(TreeTemplate<Node>& tree, const double maxLength)
{
    assert(maxLength > 0 && "Invalid maximum edge length.");

    std::vector<BeagleTreeLikelihood::AttachmentLocation> result;
    result.reserve(tree.getNumberOfNodes());

    for(bpp::Node* node : onlineAvailableEdges(tree)) {
        const double l = node->getDistanceToFather();
        // Number of points to sample on the edge - at least one
        const size_t n = static_cast<size_t>(l / maxLength) + 1;

        const double stepSize = l / (n + 1);
        for(size_t i = 0; i < n; i++) {
            result.emplace_back(node, i * stepSize);
        }
    }

    return result;
}

/// \brief sum likelihoods associated with the same node
std::unordered_map<Node*, double> accumulatePerEdgeLikelihoods(std::vector<BeagleTreeLikelihood::AttachmentLocation>& locs,
                                                               const std::vector<double>& logWeights)
{
    assert(locs.size() == logWeights.size() && "vectors differ in length");
    std::unordered_map<bpp::Node*, double> logProbByNode;
    auto locIt = locs.cbegin(), locEnd = locs.cend();
    auto logWeightIt = logWeights.cbegin();
    for(; locIt != locEnd; locIt++, logWeightIt++) {
        Node* node = locIt->first;
        auto it = logProbByNode.find(node);
        if(it == logProbByNode.end()) {
            logProbByNode[node] = *logWeightIt;
        } else {
            it->second = logSum(it->second, *logWeightIt);
        }
    }

    // Normalize
    double totalDensity = -std::numeric_limits<double>::max();
    for(auto it = logProbByNode.cbegin(), end = logProbByNode.cend(); it != end; ++it)
        totalDensity = logSum(totalDensity, it->second);
    for(auto it = logProbByNode.begin(), end = logProbByNode.end(); it != end; ++it)
        it->second -= totalDensity;
    return logProbByNode;
}

GuidedOnlineAddSequenceMove::GuidedOnlineAddSequenceMove(CompositeTreeLikelihood& calculator,
                                                         const vector<string>& taxaToAdd,
                                                         const vector<double>& proposePendantBranchLengths,
                                                         const double maxLength) :
    OnlineAddSequenceMove(calculator, taxaToAdd),
    proposePendantBranchLengths(proposePendantBranchLengths),
    maxLength(maxLength)
{
    assert(!proposePendantBranchLengths.empty() && "No proposal branch lengths!");
}


/// Choose an edge for insertion
/// This is a guided move - we calculate the likelihood with the sequence inserted at the middle of each
/// edge, then select an edge by sampling from the multinomial distribution weighted by the edge-likelihoods.
///
/// Likelihoods are calculated for each branch length in #proposePendantBranchLengths. The likelihoods from the branch
/// length with the highest / overall likelihood are used for the proposal.
///
/// When proposing by length, a collection of evenly-spaced points on the tree are created
/// This makes the proposal distribution a bit complicated - an edge may be present in the proposal set more than once.
/// accumulatePerEdgeLikelihoods (above) takes care of summing likelihoods.
const pair<Node*, double> GuidedOnlineAddSequenceMove::chooseEdge(TreeTemplate<Node>& tree, const std::string& leaf_name,
                                                                  smc::rng* rng)
{
    std::vector<BeagleTreeLikelihood::AttachmentLocation> locs = discretizeTree(tree, maxLength);
    const std::vector<std::vector<double>> attach_log_likes_by_pendant =
        calculator.calculateAttachmentLikelihoods(leaf_name, locs, proposePendantBranchLengths);
    std::vector<double> attach_log_likes(locs.size());
    auto max_double = [](const std::vector<double>& v) { return *std::max_element(v.cbegin(), v.cend()); };
    std::transform(attach_log_likes_by_pendant.cbegin(),
                   attach_log_likes_by_pendant.cend(),
                   attach_log_likes.begin(),
                   max_double);

    const std::unordered_map<bpp::Node*, double> node_log_weights = accumulatePerEdgeLikelihoods(locs, attach_log_likes);

    WeightedSelector<bpp::Node*> node_selector;
    for(auto& p : node_log_weights) {
        assert(node_log_weights.count(p.first) == 1);
        node_selector.push_back(p.first, std::exp(p.second));
    }
    assert(node_selector.size() == node_log_weights.size());

    bpp::Node* n = node_selector.choice();
    return pair<Node*,double>(n, node_log_weights.at(n));
}

/// Propose branch-lengths around ML value
TripodOptimizer GuidedOnlineAddSequenceMove::optimizeBranchLengths(const Node* insertEdge,
                                                                   const std::string& newLeafName,
                                                                   double& distalBranchLength,
                                                                   double& pendantBranchLength)
{
    const double d = insertEdge->getDistanceToFather();

    TripodOptimizer optim = calculator.createOptimizer(insertEdge, newLeafName);

    double pendant = 1e-8;
    double distal = d / 2;

    // Optimize distal, pendant up to 5 times
    for(size_t i = 0; i < 5; i++) {
        const double newDistal = optim.optimizeDistal(distal, pendant);
        if(std::abs(newDistal - distal) < TripodOptimizer::TOLERANCE)
            break;

        const double newPendant = optim.optimizePendant(newDistal, pendant);
        if(std::abs(newPendant - pendant) < TripodOptimizer::TOLERANCE)
            break;

        pendant = newPendant;
        distal = newDistal;
    }

    distalBranchLength = distal;
    pendantBranchLength = pendant;

    return optim;
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
    return AttachmentProposal { n, edgeLogDensity, distal, distalLogDensity, pendantBranchLength, pendantLogDensity, mlDistal, mlPendant, false };
}

}} // namespaces

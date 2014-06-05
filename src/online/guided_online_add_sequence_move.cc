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
std::vector<BeagleTreeLikelihood::AttachmentLocation> discretizeTree(TreeTemplate<Node>& tree, const size_t n)
{
    // First we need the total tree length
    std::vector<bpp::Node*> nodes = onlineAvailableEdges(tree);
    const auto f = [](const double acc, const bpp::Node* n) { return acc + n->getDistanceToFather(); };
    const double totalTreeLength = std::accumulate(nodes.begin(), nodes.end(), 0.0, f);

    assert(n > 1 && "Invalid number of partitions");
    const double stepSize = totalTreeLength / (static_cast<double>(n) + 1);

    auto nodeIt = nodes.begin(), nodeEnd = nodes.end();

    std::vector<BeagleTreeLikelihood::AttachmentLocation> result;
    result.reserve(n);
    // the amount of total branch length covered up to the node pointed to by nodeIt.
    double lengthCoveredToNode = 0;

    for(size_t i = 1; i <= n; i++) {
        // Target position along the tree's total branch length
        const double target = stepSize * i;

        // Advance nodeIt until the range:
        // [lengthCoveredToNode, lengthCoveredToNode + (*nodeIt)->getDistanceToFather()]
        // contains target.
        while(lengthCoveredToNode + (*nodeIt)->getDistanceToFather() < target) {
            lengthCoveredToNode += (*nodeIt)->getDistanceToFather();
            ++nodeIt;
            assert(nodeIt != nodeEnd && "ran out of nodes");
            if(!(*nodeIt)->hasDistanceToFather())
                ++nodeIt;
            assert(nodeIt != nodeEnd && "ran out of nodes");
        }

        // We've advanced to the correct node
        result.emplace_back(*nodeIt, target - lengthCoveredToNode);
        assert(result.back().second <= (*nodeIt)->getDistanceToFather());
    }

    assert(result.size() == n && "Result size does not match expected");

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
    double totalDensity = 0;
    for(auto it = logProbByNode.cbegin(), end = logProbByNode.cend(); it != end; ++it)
        totalDensity = logSum(totalDensity, it->second);
    for(auto it = logProbByNode.begin(), end = logProbByNode.end(); it != end; ++it)
        it->second -= totalDensity;
    return logProbByNode;
}

GuidedOnlineAddSequenceMove::GuidedOnlineAddSequenceMove(CompositeTreeLikelihood& calculator,
                                                         const vector<string>& taxaToAdd,
                                                         const vector<double>& proposePendantBranchLengths,
                                                         const bool byLength) :
    OnlineAddSequenceMove(calculator, taxaToAdd),
    proposePendantBranchLengths(proposePendantBranchLengths),
    byLength(byLength)
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
    if(byLength) {
        std::vector<BeagleTreeLikelihood::AttachmentLocation> locs = discretizeTree(tree, tree.getNumberOfNodes() * 2);
        const std::vector<std::vector<double>> attach_log_likes_by_pendant =
            calculator.calculateAttachmentLikelihoods(leaf_name, locs, proposePendantBranchLengths);
        std::vector<double> attach_log_likes(locs.size());
        auto max_double = [](const std::vector<double>& v) { return *std::max_element(v.begin(), v.end()); };
        std::transform(attach_log_likes_by_pendant.cbegin(),
                       attach_log_likes_by_pendant.cend(),
                       attach_log_likes.begin(),
                       max_double);

        std::unordered_map<bpp::Node*, double> node_log_weights = accumulatePerEdgeLikelihoods(locs, attach_log_likes);

        WeightedSelector<bpp::Node*> node_selector;
        for(auto& p : node_log_weights) {
            node_selector.push_back(p.first, p.second);
        }

        bpp::Node* n = node_selector.choice();
        return pair<Node*,double>(n, node_log_weights.at(n));
    } else {
        // By edge
        const std::vector<double> edge_log_likes = calculator.edgeLogLikelihoods(leaf_name, proposePendantBranchLengths);

        const double totalLike = std::accumulate(edge_log_likes.begin() + 1, edge_log_likes.end(),
                                                 edge_log_likes[0],
                                                 [](double a, double b) { return logSum(a, b); });
        vector<double> edge_likes(edge_log_likes.size());
        std::transform(edge_log_likes.begin(), edge_log_likes.end(), edge_likes.begin(),
                       [&totalLike](double d) { return std::exp(d - totalLike); });

        // Select an edge
        std::vector<unsigned> indexes(edge_likes.size());
        rng->Multinomial(1, edge_likes.size(), edge_likes.data(), indexes.data());
        // Only one value should be selected - find it
        auto positive = [](const unsigned x) { return x > 0; };
        std::vector<unsigned>::const_iterator it = std::find_if(indexes.begin(), indexes.end(),
                                                                positive);
        assert(it != indexes.end());
        const size_t idx = it - indexes.begin();

        Node* n = onlineAvailableEdges(tree)[idx];
        return pair<Node*,double>(n, edge_log_likes[idx] - totalLike);
    }
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

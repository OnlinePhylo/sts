#include "guided_online_add_sequence_move.h"
#include "tree_particle.h"
#include "beagle_tree_likelihood.h"
#include "composite_tree_likelihood.h"
#include "likelihood_vector.h"
#include "tripod_optimizer.h"
#include "util.h"
#include "online_util.h"

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

GuidedOnlineAddSequenceMove::GuidedOnlineAddSequenceMove(CompositeTreeLikelihood& calculator,
                                                         const vector<string>& taxaToAdd,
                                                         const vector<double>& proposePendantBranchLengths) :
    OnlineAddSequenceMove(calculator, taxaToAdd),
    proposePendantBranchLengths(proposePendantBranchLengths)
{
    assert(!proposePendantBranchLengths.empty() && "No proposal branch lengths!");
}

/// \brief Discretize the tree into \c n evenly spaced points
std::vector<BeagleTreeLikelihood::AttachmentLocation> discretizeTree(TreeTemplate<Node>& tree, const size_t n)
{
    // First we need the total tree length
    std::vector<bpp::Node*> nodes = onlineAvailableEdges(tree);
    auto f = [](const double acc, const bpp::Node* n) { return acc + n->getDistanceToFather(); };
    const double totalTreeLength = std::accumulate(nodes.begin(), nodes.end(), 0.0, f);

    for(const bpp::Node* n : nodes) {
        assert(n != nullptr && "Null node?");
    }

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

/// Choose an edge for insertion
/// This is a guided move - we calculate the likelihood with the sequence inserted at the middle of each
/// edge, then select an edge by sampling from the multinomial distribution weighted by the edge-likelihoods.
///
/// Likelihoods are calculated for each branch length in #proposePendantBranchLengths. The likelihoods from the branch
/// length with the highest / overall likelihood are used for the proposal.
const pair<Node*, double> GuidedOnlineAddSequenceMove::chooseEdge(TreeTemplate<Node>& tree, const std::string& leaf_name,
                                                                  smc::rng* rng)
{
    // Alright alright - try discretizing the tree.
    // TODO: hard-coded constant is temporary.
    const std::vector<BeagleTreeLikelihood::AttachmentLocation> locs = discretizeTree(tree, 100);
    const std::vector<double> attach_log_likes = calculator.calculateAttachmentLikelihoods(leaf_name, locs);

    const double bestLogLike = *std::max_element(attach_log_likes.begin(), attach_log_likes.end());
    vector<double> attach_likes(attach_log_likes.size());
    std::transform(attach_log_likes.begin(), attach_log_likes.end(), attach_likes.begin(),
                   [&bestLogLike](double p) { return std::exp(p - bestLogLike); });

    // Select an edge
    std::vector<unsigned> indexes(attach_likes.size());
    rng->Multinomial(1, attach_likes.size(), attach_likes.data(), indexes.data());
    // Only one value should be selected - find it
    auto positive = [](const unsigned x) { return x > 0; };
    std::vector<unsigned>::const_iterator it = std::find_if(indexes.begin(), indexes.end(),
                                                            positive);
    assert(it != indexes.end());
    const size_t idx = it - indexes.begin();

    Node* n = locs.at(idx).first;
    return pair<Node*,double>(n, attach_log_likes[idx] - bestLogLike);
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

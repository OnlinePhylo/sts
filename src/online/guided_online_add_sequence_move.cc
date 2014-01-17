#include "guided_online_add_sequence_move.h"
#include "tree_particle.h"
#include "beagle_tree_likelihood.h"
#include "composite_tree_likelihood.h"
#include "likelihood_vector.h"
#include "tripod_optimizer.h"
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
const pair<Node*, double> GuidedOnlineAddSequenceMove::chooseEdge(TreeTemplate<Node>& tree, const std::string& leaf_name,
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
TripodOptimizer GuidedOnlineAddSequenceMove::optimizeBranchLengths(const Node* insertEdge,
                                                                   const std::string& newLeafName,
                                                                   double& distalBranchLength,
                                                                   double& pendantBranchLength)
{
    const double d = insertEdge->getDistanceToFather();

    std::shared_ptr<BeagleTreeLikelihood> btl(calculator.calculator());
    assert(btl->freeBufferCount() >= 2);

    TripodOptimizer optim(btl, insertEdge, newLeafName, d);

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

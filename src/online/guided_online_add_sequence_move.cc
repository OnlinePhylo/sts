#include "guided_online_add_sequence_move.h"
#include "tree_particle.h"
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

#include <gsl/gsl_cdf.h>

using namespace std;
using namespace bpp;

namespace sts { namespace online {

std::vector<AttachmentLocation> divideEdge(Node* node,
                                           const double maxLength)
{
    assert(node != nullptr);
    assert(node->hasDistanceToFather());
    std::vector<AttachmentLocation> result;
    const double l = node->getDistanceToFather();
    // Number of points to sample on the edge - at least one
    const size_t n = static_cast<size_t>(l / maxLength) + 1;

    const double stepSize = l / (n + 1);
    for(size_t i = 1; i <= n; i++) {
        result.emplace_back(node, i * stepSize);
    }
    return result;
}

/// \brief Generate attachment locations for `nodes`
///
/// Divides the edges in the input node list such that:
///
/// * There is at least one attachment attempt on each edge (by default: at the midpoint)
/// * the distance from any node to an attachment location in both directions is always `<= maxLength`
std::vector<AttachmentLocation> divideEdges(std::vector<Node*>& nodes,
                                            const double maxLength)
{
    assert(maxLength > 0.0 && "Invalid maximum edge length");
    std::vector<AttachmentLocation> result;
    result.reserve(nodes.size());

    for(bpp::Node* node : nodes) {
        std::vector<AttachmentLocation> innerResult = divideEdge(node, maxLength);
        std::copy(innerResult.begin(), innerResult.end(), std::back_inserter(result));
    }

    return result;
}

/// \brief Discretize the tree into \c n evenly spaced points
std::vector<AttachmentLocation> divideTreeEdges(TreeTemplate<Node>& tree, const double maxLength)
{
    std::vector<Node*> nodes = onlineAvailableEdges(tree);
    return divideEdges(nodes, maxLength);
}


/// \brief average likelihoods associated with the same node
///
/// \param locs Attachment locations. Nodes may be repeated, in which case their likelihoods are averaged.
/// \param logWeights log-likelihood associated with each location in `locs`.
/// \return A normalized probability attaching to each edge.
/// \pre `locs.size() == logWeights.size()`
std::vector<std::pair<bpp::Node*, double> > GuidedOnlineAddSequenceMove::accumulatePerEdgeLikelihoods(std::vector<AttachmentLocation>& locs,
                                                               const std::vector<double>& logWeights) const
{
    assert(locs.size() == logWeights.size() && "vectors differ in length");
    std::vector<std::pair<bpp::Node*, double> > logProbByNode;
    std::unordered_map<bpp::Node*, int> countByNode;
    auto locIt = locs.cbegin(), locEnd = locs.cend();
    auto logWeightIt = logWeights.cbegin();

    // Sum likelihood by node
    for(; locIt != locEnd; locIt++, logWeightIt++) {
        Node* node = locIt->node;
        auto it = std::find_if(logProbByNode.begin(), logProbByNode.end(),
                               [node](const std::pair<bpp::Node*, double>& p){return p.first == node;});
        if(it == logProbByNode.end()) {
            logProbByNode.push_back(std::make_pair(node, *logWeightIt));
            countByNode[node] = 1;
        } else {
            size_t pos = it - logProbByNode.begin();
            logProbByNode[pos].second = logSum(logProbByNode[pos].second, *logWeightIt);
            countByNode[node] += 1;
        }
    }

    // Calculate mean likelihood per node
    for(auto &p : logProbByNode) {
        // Number of sampled points
        const int count = countByNode.at(p.first);
        assert(count > 0);
        if(count > 1)
            p.second -= std::log(static_cast<double>(count));
    }

    // Total likelihood
    double totalDensity = -std::numeric_limits<double>::max();
    for(auto it = logProbByNode.begin(), end = logProbByNode.end(); it != end; ++it){
    	it->second *= _heating;
        totalDensity = logSum(totalDensity, it->second);
    }

    // Normalize node likelihoods
    for(auto it = logProbByNode.begin(), end = logProbByNode.end(); it != end; ++it)
        it->second -= totalDensity;
    return logProbByNode;
}

GuidedOnlineAddSequenceMove::GuidedOnlineAddSequenceMove(CompositeTreeLikelihood& calculator,
                                                         const std::vector<std::string>& sequenceNames,
                                                         const vector<string>& taxaToAdd,
                                                         const vector<double>& proposePendantBranchLengths,
                                                         const double maxLength,
                                                         const size_t subdivideTop) :
    OnlineAddSequenceMove(calculator, sequenceNames, taxaToAdd),
    proposePendantBranchLengths(proposePendantBranchLengths),
    maxLength(maxLength),
    subdivideTop(subdivideTop)
{
    _proposalMethodName = "GuidedOnlineAddSequenceMove";
    assert(!proposePendantBranchLengths.empty() && "No proposal branch lengths!");
    _heating = 0.05;
}

/// Passing by value purposefully here
std::vector<std::pair<bpp::Node*, double> > GuidedOnlineAddSequenceMove::subdivideTopN(std::vector<AttachmentLocation> locs,
                                                                             const std::vector<double>& logWeights,
                                                                             const std::string& leafName)
{
    assert(logWeights.size() == locs.size() && "Invalid size");

    const size_t subdivideTop = std::min(this->subdivideTop, logWeights.size());

    // Lookup table for log weights
    std::unordered_map<const bpp::Node*, double> nodeLogWeights;
    nodeLogWeights.reserve(locs.size());
    {
        auto locIt = locs.begin(), locEnd = locs.end();
        auto weightIt = logWeights.begin();
        for(; locIt != locEnd; ++locIt, ++weightIt) {
            assert(nodeLogWeights.find(locIt->node) == nodeLogWeights.end());
            nodeLogWeights[locIt->node] = *weightIt;
        }
    }

    // Sort locations by *descending* likelihood
    auto key = [&nodeLogWeights](const AttachmentLocation& x,
                                 const AttachmentLocation& y) -> bool {
        return nodeLogWeights.at(x.node) > nodeLogWeights.at(y.node);
    };
    std::sort(locs.begin(), locs.end(), key);

    std::vector<AttachmentLocation> tmpLocs;
    std::vector<double> tmpLogLikes;
    for(size_t i = 0; i < locs.size(); i++) {
        AttachmentLocation& loc = locs.at(i);
        if(i >= subdivideTop || loc.node->getDistanceToFather() < 2 * maxLength) {
            // Edge does not need to be divided further
            tmpLocs.push_back(loc);
            tmpLogLikes.push_back(nodeLogWeights.at(loc.node));
        } else {
            // Edge needs to be divided further
            assert(loc.node != nullptr);
            for(AttachmentLocation& l : divideEdge(loc.node, maxLength)) {
                tmpLocs.push_back(l);
                // Posterior
                std::vector<double> ll;
                for(double pendantLength : proposePendantBranchLengths){
                    ll.push_back(calculator(*l.node, leafName, pendantLength, l.distal, l.node->getDistanceToFather()-l.distal));
                }
                
                tmpLogLikes.push_back(*std::max_element(ll.begin(), ll.end()));
            }
        }
    }
    assert(tmpLocs.size() >= locs.size());

    // Update
    return accumulatePerEdgeLikelihoods(tmpLocs, tmpLogLikes);
}

/// Choose an edge for insertion
/// This is a guided move - we calculate the likelihood with the sequence inserted at the middle of each
/// edge, then select an edge by sampling from the multinomial distribution weighted by the edge-likelihoods.
///
/// Likelihoods are calculated for each branch length in #proposePendantBranchLengths. The likelihoods from the branch
/// length with the highest / overall likelihood are used for the proposal.
///
/// When proposing by length, edges are divided such that unsampled segments are no longer than `maxLength`
/// This makes the proposal distribution a bit complicated - an edge may be present in the proposal set more than once.
/// accumulatePerEdgeLikelihoods (above) takes care of averaging likelihoods.
const pair<Node*, double> GuidedOnlineAddSequenceMove::chooseEdge(TreeTemplate<Node>& tree,
                                                                  const std::string& leafName,
                                                                  smc::rng* rng, size_t particleID)
{
    // If subdivideTop is set, we do not subdivide edges here, rather
    // divide in half once and subdivide the top N edges later
    std::vector<AttachmentLocation> locs =
        divideTreeEdges(tree,
                        (subdivideTop > 0.0 ?
                         std::numeric_limits<double>::max() :
                         maxLength));
    
    std::vector<std::vector<double>> attachLogLikesByPendant;
    for(AttachmentLocation& l : locs) {
        // Posterior
        std::vector<double> ll;
        for(double pendantLength : proposePendantBranchLengths){
            ll.push_back(calculator(*l.node, leafName, pendantLength, l.distal, l.node->getDistanceToFather()-l.distal));
        }
        
        attachLogLikesByPendant.push_back(ll);
    }

    std::vector<double> attachLogLikes(locs.size());

    auto maxDouble = [](const std::vector<double>& v) {
        assert(!v.empty());
        return *std::max_element(v.cbegin(), v.cend());
    };

    std::transform(attachLogLikesByPendant.cbegin(),
                   attachLogLikesByPendant.cend(),
                   attachLogLikes.begin(),
                   maxDouble);

    std::vector<std::pair<bpp::Node*, double> > nodeLogWeights;

    if(subdivideTop > 0) {
        // Hybrid scheme
        nodeLogWeights = subdivideTopN(locs, attachLogLikes, leafName);
    } else {
        nodeLogWeights = accumulatePerEdgeLikelihoods(locs, attachLogLikes);
    }
    std::vector<std::pair<size_t, double>> probabilities;
    probabilities.reserve(nodeLogWeights.size());
    
    WeightedSelector<bpp::Node*> nodeSelector{*rng};
    for(auto& p : nodeLogWeights) {
        double prob = exp(p.second);
        nodeSelector.push_back(p.first, prob);
        probabilities.push_back(make_pair(p.first->getId(), prob));
    }
    assert(nodeSelector.size() == nodeLogWeights.size());
    
    _probs[particleID] = probabilities;

    bpp::Node* n = nodeSelector.choice();
    auto it = std::find_if( nodeLogWeights.begin(), nodeLogWeights.end(),
                           [n](const std::pair<bpp::Node*, double>& element){return element.first == n;});
    return pair<Node*,double>(n, it->second);
}

/// Propose branch-lengths around ML value
void GuidedOnlineAddSequenceMove::optimizeBranchLengths(const Node* insertEdge,
                                                                   const std::string& newLeafName,
                                                                   double& distalBranchLength,
                                                                   double& pendantBranchLength)
{
    const double d = insertEdge->getDistanceToFather();
    TripodOptimizer optim(calculator, insertEdge, newLeafName, d);

    double pendant = 1e-8;
    double distal = d / 2;

    // Optimize distal, pendant up to 5 times
    for(size_t i = 0; i < 5; i++) {
        size_t nChanged = 0;

        // Optimize distal if it's longer than the tolerance
        const double newDistal = (d <= TripodOptimizer::TOLERANCE) ?
            distal :
            optim.optimizeDistal(distal, pendant);

        if(std::abs(newDistal - distal) > TripodOptimizer::TOLERANCE)
            nChanged++;

        const double newPendant = optim.optimizePendant(newDistal, pendant);
        if(std::abs(newPendant - pendant) > TripodOptimizer::TOLERANCE)
            nChanged++;

        pendant = newPendant;
        distal = newDistal;

        if(!nChanged)
            break;
    }

    distalBranchLength = distal;
    pendantBranchLength = pendant;
}

/// Distal branch length proposal
/// We propose from Gaussian(mlDistal, edgeLength / 4) with support truncated to [0, edgeLength]
std::pair<double, double> GuidedOnlineAddSequenceMove::proposeDistal(bpp::Node& n, const std::string& leafName, const double mlDistal, const double mlPendant, smc::rng* rng) const
{
    assert(mlDistal <= n.getDistanceToFather());
    const double edgeLength = n.getDistanceToFather();
    
    double distal = -1;
    // HACK/ARBITRARY: proposal standard deviation
    const double sigma = edgeLength / 4;

    // Handle very small branch lengths - attach with distal BL of 0
    if(edgeLength < 1e-8)
        distal = 0;
    else {
        do {
            distal = rng->NormalTruncated(mlDistal, sigma, 0.0);
        } while(distal < 0 || distal > edgeLength);
    }
    assert(!std::isnan(distal));

    // Log density: for CDF F(x) and PDF g(x), limited to the interval (a, b]:
    //
    // g'(x) =   g(x) / [F(b) - F(a)]
    //
    // We are limited to (0, d].
    //
    // GSL gaussian CDFs are for mean 0, hence the mlDistal substraction here.
    const double distalLogDensity = std::log(gsl_ran_gaussian_pdf(distal - mlDistal, sigma)) -
        std::log(gsl_cdf_gaussian_P(edgeLength - mlDistal, sigma) - gsl_cdf_gaussian_P(- mlDistal, sigma));
    assert(!std::isnan(distalLogDensity));

    return std::pair<double, double>(distal, distalLogDensity);
}

std::pair<double, double> GuidedOnlineAddSequenceMove::proposePendant(bpp::Node& n, const std::string& leafName, const double mlPendant, const double distalBranchLength, smc::rng* rng) const
{
    const double pendantBranchLength = rng->Exponential(mlPendant);
    const double pendantLogDensity = std::log(gsl_ran_exponential_pdf(pendantBranchLength, mlPendant));
    return std::pair<double, double>(pendantBranchLength, pendantLogDensity);
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
    
    // Step 1: Select attachment branch
    Node* n = nullptr;
    double edgeLogDensity;
    size_t toAddCount = std::distance(taxaToAdd.begin(),taxaToAdd.end());
    if(_toAddCount == toAddCount && _probs.find(value->particleID) != _probs.end()){
        std::vector<bpp::Node*> nodes = onlineAvailableEdges(*tree);

        const std::vector<std::pair<size_t, double>>& probabilities = _probs[value->particleID];
        WeightedSelector<bpp::Node*> selector{*rng};
        for(auto& pair: probabilities){
            size_t idx = pair.first;
            auto it = std::find_if( nodes.begin(), nodes.end(),
                                   [idx](const bpp::Node* element){return element->getId() == idx;});
            selector.push_back(*it, pair.second);
        }
        n = selector.choice();
        int idx = n->getId();
        auto it = std::find_if( probabilities.begin(), probabilities.end(),
                               [idx](const std::pair<size_t, double>& element){return element.first == idx;});
        edgeLogDensity = log(it->second);
    }
    else{
        std::tie(n, edgeLogDensity) = chooseEdge(*tree, leafName, rng, value->particleID);
    }
    
//    calculator.calculateAttachmentLikelihood(leafName, n, 0);
    
    // Calculate MLEs of distal and pendant branch lengths
    bool found = false;
    std::unordered_map<size_t, std::unordered_map<size_t, std::pair<double,double>>>::const_iterator iter= _mles.find(value->particleID);
    if(_toAddCount == toAddCount && iter != _mles.cend()){
        std::unordered_map<size_t, std::pair<double,double>>::const_iterator iter2 = iter->second.find(n->getId());
        if(iter2 != iter->second.cend()){
            std::tie(_mleDistal, _mlePendant) = iter2->second;
            found = true;
        }
    }
    if(!found){
        optimizeBranchLengths(n, leafName, _mleDistal, _mlePendant);
        _mles[value->particleID][n->getId()] = std::make_pair(_mleDistal, _mlePendant);
    }


    // Step2:  proposal distal branch length
    double distalBranchLength, distalLogDensity;
    std::tie(distalBranchLength, distalLogDensity) = proposeDistal(*n, leafName, _mleDistal, _mlePendant, rng);

    // Step 3: propose pendant branch length
    double pendantBranchLength, pendantLogDensity;
    std::tie(pendantBranchLength, pendantLogDensity) = proposePendant(*n, leafName, _mlePendant, distalBranchLength, rng);
    assert(!std::isnan(pendantLogDensity));
    
    return AttachmentProposal { n, edgeLogDensity, distalBranchLength, distalLogDensity, pendantBranchLength, pendantLogDensity, _mleDistal, _mlePendant, _proposalMethodName };
}

}} // namespaces

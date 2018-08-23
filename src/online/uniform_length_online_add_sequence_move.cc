#include "uniform_length_online_add_sequence_move.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <memory>
#include <numeric>

#include "composite_tree_likelihood.h"
#include "online_util.h"
#include "tree_particle.h"
#include "util.h"
#include "weighted_selector.h"

using namespace std;
using namespace bpp;

namespace sts { namespace online {

// called when --proposal cmdline arg is "uniform-length"
//
UniformLengthOnlineAddSequenceMove::UniformLengthOnlineAddSequenceMove(std::vector<std::unique_ptr<CompositeTreeLikelihood>>& calculator,
                                                                       const std::vector<std::string>& sequenceNames,
                                                                       const vector<string>& taxaToAdd,
                                                                       std::function<std::pair<double,double>(smc::rng*)> branchLengthProposer) :
    OnlineAddSequenceMove(calculator, sequenceNames, taxaToAdd),
    branchLengthProposer(branchLengthProposer)
{ }


// proposes a distal position and pendant length.
// chooses an edge from a multinomial distribution weighted by length of the edges
// distal position is selected from a uniform distrbution across the edge length.
// pendant length is extracted from branchLengthProposer passed in at instace creation.
// 
AttachmentProposal UniformLengthOnlineAddSequenceMove::propose(const std::string&, smc::particle<TreeParticle>& particle, smc::rng* rng)
{
    TreeParticle* value = particle.GetValuePointer();
    unique_ptr<TreeTemplate<bpp::Node>>& tree = value->tree;

    std::vector<bpp::Node*> nodes = onlineAvailableEdges(*tree);
    std::vector<double> lengths;
    lengths.reserve(nodes.size());

    // select an edge from a multinomial distribution weighted by length of the edges
    WeightedSelector<size_t> selector{*rng};
    size_t i = 0;
    for(bpp::Node* n : nodes) {
	// csw - remove some logic that guarded against root and rightOfRoot.
	// onlineAvailableEdges() already filters those out.
	lengths.push_back(n->getDistanceToFather());
	selector.push_back(i++, n->getDistanceToFather());
    }
    assert(selector.size() > 0 && "No eligible nodes!");

    i = selector.choice();
    Node* n = nodes.at(i);

    // branch lengths
    double pendant, pendantLogDensity;
    std::tie(pendant, pendantLogDensity) = branchLengthProposer(rng);
    const double d = n->getDistanceToFather();
    const double distal = rng->UniformS() * d;

    // Prior probability of the edge proposal is the relative branch length
    const double totalLength = std::accumulate(lengths.begin(), lengths.end(), 0.0);
    const double edgePriorDensity = std::log(lengths[i]) - std::log(totalLength);
    return AttachmentProposal {n, edgePriorDensity, distal, 0.0, pendant, pendantLogDensity, 0.0, 0.0, "UniformLengthOnlineAddSequenceMove"};
}

}} // namespaces

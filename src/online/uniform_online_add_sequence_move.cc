#include "uniform_online_add_sequence_move.h"
#include "tree_particle.h"
#include "beagle_tree_likelihood.h"
#include "composite_tree_likelihood.h"
#include "util.h"

#include <algorithm>
#include <iterator>
#include <cassert>
#include <cmath>
#include <memory>

using namespace std;
using namespace bpp;
using sts::util::beagle_check;

namespace sts { namespace online {

UniformOnlineAddSequenceMove::UniformOnlineAddSequenceMove(CompositeTreeLikelihood& calculator,
                                                           const std::vector<std::string>& sequenceNames,
                                                           const vector<string>& taxaToAdd,
                                                           std::function<std::pair<double,double>(smc::rng*)> branchLengthProposer) :
    OnlineAddSequenceMove(calculator, sequenceNames, taxaToAdd),
    branchLengthProposer(branchLengthProposer)
{ }

AttachmentProposal UniformOnlineAddSequenceMove::propose(const std::string& leafName, smc::particle<TreeParticle>& particle, smc::rng* rng)
{
    TreeParticle* value = particle.GetValuePointer();
    unique_ptr<TreeTemplate<bpp::Node>>& tree = value->tree;

    size_t idx = rng->UniformDiscrete(0, tree->getNumberOfNodes() - 3);
    std::vector<bpp::Node*> nodes = tree->getNodes();
    while(nodes[idx] == tree->getRootNode()->getSon(1) || nodes[idx] == tree->getRootNode())
        idx++;

    Node* n = nodes[idx];

    // branch lengths
    double pendant, pendantLogDensity;
    std::tie(pendant, pendantLogDensity) = branchLengthProposer(rng);
    const double d = n->getDistanceToFather();
    const double distal = rng->UniformS() * d;

    return AttachmentProposal {n, 0.0, distal, 0.0, pendant, pendantLogDensity, 0.0, 0.0, "UniformOnlineAddSequenceMove"};
}

}} // namespaces

#include "online_add_sequence_move.h"
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

OnlineAddSequenceMove::OnlineAddSequenceMove(CompositeTreeLikelihood& calculator,
                                             const vector<string>& taxaToAdd) :
    calculator(calculator),
    taxaToAdd(std::begin(taxaToAdd), std::end(taxaToAdd)),
    lastTime(-1)
{ }

void OnlineAddSequenceMove::addProposalRecord(const ProposalRecord& proposalRecord)
{
    proposalRecords_.push_back(proposalRecord);
}

const std::vector<ProposalRecord> OnlineAddSequenceMove::getProposalRecords() const
{
    return proposalRecords_;
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

    AttachmentProposal proposal = propose(taxaToAdd.front(), particle, rng);

    // New internal node, new leaf
    Node* new_node = new Node(tree->getNumberOfNodes());
    Node* new_leaf = new Node(new_node->getId() + 1, taxaToAdd.front());
    new_node->addSon(new_leaf);
    new_leaf->setDistanceToFather(proposal.pendantBranchLength);

    assert(proposal.edge->hasFather());
    Node* father = proposal.edge->getFather();

    // branch lengths
    const double d = proposal.edge->getDistanceToFather();
    assert(proposal.distalBranchLength <= d && "Distal branch length exceeds total!");

    // Swap `new_node` in for `n`
    // Note: use {add,remove}Son, rather than {remove,set}Father -
    // latter functions do not update parent sons list.
    father->addSon(new_node);
    father->removeSon(proposal.edge);
    new_node->addSon(proposal.edge);

    // Attachment branch lengths
    new_node->setDistanceToFather(d - proposal.distalBranchLength);
    proposal.edge->setDistanceToFather(proposal.distalBranchLength);

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

    const double log_like = calculator();

    const double orig_weight = particle.GetLogWeight();
    particle.AddToLogWeight(log_like);
    particle.AddToLogWeight(-proposal.logProposalDensity());
    particle.AddToLogWeight(-orig_ll);
    const double new_weight = particle.GetLogWeight();

    assert(!std::isnan(particle.GetLogWeight()));

    addProposalRecord({time, orig_ll, log_like, orig_weight, new_weight, proposal});
}

}} // namespaces

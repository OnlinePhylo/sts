#include "online_add_sequence_move.h"
#include "tree_particle.h"
#include "composite_tree_likelihood.h"
#include "util.h"

#include <algorithm>
#include <iterator>
#include <cassert>
#include <cmath>
#include <memory>

#include <gsl/gsl_randist.h>

using namespace std;
using namespace bpp;

namespace sts { namespace online {

OnlineAddSequenceMove::OnlineAddSequenceMove(CompositeTreeLikelihood& calculator,
                                             const vector<string>& sequenceNames,
                                             const vector<string>& taxaToAdd) :
    calculator(calculator),
    _sequenceNames(sequenceNames),
    taxaToAdd(std::begin(taxaToAdd), std::end(taxaToAdd)),
    _toAddCount(-1),
    _counter(0),
    lastTime(-1),
    _empiricalGammaProposal(nullptr),
    _empiricalMVNProposal(nullptr)
{ }
    
void OnlineAddSequenceMove::setGammaProposal(std::function<std::tuple<double, double>(smc::rng*)> empiricalGammaProposal){
    _empiricalGammaProposal = empiricalGammaProposal;
}

void OnlineAddSequenceMove::setMVNProposal(std::function<std::tuple<std::map<std::string, double>, double>(smc::rng*)> empiricalMVNProposal){
    _empiricalMVNProposal = empiricalMVNProposal;
}

void OnlineAddSequenceMove::addProposalRecord(const ProposalRecord& proposalRecord)
{
    proposalRecords_.push_back(proposalRecord);
}

const std::vector<ProposalRecord> OnlineAddSequenceMove::getProposalRecords() const
{
    return proposalRecords_;
}
    
void OnlineAddSequenceMove::addTaxa(const vector<string>& taxa){
    _toAddCount = -1;
    lastTime = -1;
    taxaToAdd.clear();
    auto it = taxaToAdd.begin();
    taxaToAdd.insert_after(it, taxa.cbegin(), taxa.cend());
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
    
    size_t toAddCount = std::distance(taxaToAdd.begin(),taxaToAdd.end());
    
    // we start a new iteration
    if(_toAddCount != toAddCount){
        _probs.clear();
        _mles.clear();
        _counter = 0;
    }

    AttachmentProposal proposal = propose(taxaToAdd.front(), particle, rng);
    
//    const double log_like = calculator(*proposal.edge, taxaToAdd.front(), proposal.pendantBranchLength, proposal.distalBranchLength, proposal.edge->getDistanceToFather()-proposal.distalBranchLength);
//    log_like += calculator.sumAdditionalLogLikes();
    _toAddCount = toAddCount;
    value->particleID = _counter++;

    const int new_node_id = 2*_sequenceNames.size()-1 - std::distance(taxaToAdd.cbegin(), taxaToAdd.cend());

    // New internal node, new leaf
    Node* new_node = new Node(new_node_id, "node"+std::to_string(tree->getNumberOfNodes()));
    
    size_t idx = find(_sequenceNames.begin(), _sequenceNames.end(), taxaToAdd.front()) - _sequenceNames.begin();
    assert(idx<_sequenceNames.size());
    Node* new_leaf = new Node(idx, taxaToAdd.front());
    new_node->addSon(new_leaf);
    new_leaf->setDistanceToFather(proposal.pendantBranchLength);

    assert(proposal.edge->hasFather());
    Node* father = proposal.edge->getFather();

    // branch lengths
    const double d = proposal.edge->getDistanceToFather();
    assert(proposal.distalBranchLength <= d && "Distal branch length exceeds total!");

    // Swap `new_node` in for `n`
    // Son of the root with branch length = 0 an index = 1 is not moved to the other side
    size_t pos = father->getSonPosition(proposal.edge);
    father->setSon(pos, new_node);
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
    
    double alpha;
    double alphaLogP = 0;
    if(value->rateDist->getNumberOfCategories() > 1){
        std::tie(alpha, alphaLogP) = _empiricalGammaProposal(rng);
        value->rateDist->setParameterValue("alpha", alpha);
    }
    std::map<std::string, double> parameters;
    double modelLogP = 0;
    // GTR
    if(value->model->getParameters().size() > 1){
        std::tie(parameters, modelLogP) = _empiricalMVNProposal(rng);
        for(const auto& val : parameters){
            value->model->setParameterValue(val.first, val.second);
        }
    }

    // Calculate new LL - need to re-initialize since nodes have been added
    // TODO: Should nodes be allocated dynamically?
    calculator.initialize(*value->model, *value->rateDist, *value->tree);

    const double log_like = calculator();
    value->logP = log_like;
    std::cout << log_like << std::endl;
    

    const double orig_weight = particle.GetLogWeight();
    particle.AddToLogWeight(log_like);
    particle.AddToLogWeight(-proposal.logProposalDensity() - alphaLogP - modelLogP);
    particle.AddToLogWeight(-orig_ll);
    const double new_weight = particle.GetLogWeight();

    assert(!std::isnan(particle.GetLogWeight()));

    addProposalRecord({time, orig_ll, log_like, orig_weight, new_weight, proposal});
}

}} // namespaces

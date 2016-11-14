#include "lcfit_online_add_sequence_move.h"

#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <gsl/gsl_cdf.h>
#include <lcfit_rejection_sampler.h>
#include <lcfit_select.h>

#include "attachment_likelihood.h"
#include "composite_tree_likelihood.h"
#include "guided_online_add_sequence_move.h"
#include "lcfit_rejection_sampler.h"
#include "tree_particle.h"
#include "online_util.h"
#include "weighted_selector.h"

namespace sts { namespace online {

LcfitOnlineAddSequenceMove::LcfitOnlineAddSequenceMove(CompositeTreeLikelihood& calculator,
                                                       const std::vector<std::string>& taxaToAdd,
                                                       const std::vector<double>& proposePendantBranchLengths,
                                                       const double maxLength,
                                                       const size_t subdivideTop,
                                                       const double expPriorMean) :
    GuidedOnlineAddSequenceMove(calculator, taxaToAdd, proposePendantBranchLengths, maxLength, subdivideTop),
    expPriorMean_(expPriorMean),
    lcfit_failures_(0),
    lcfit_attempts_(0)
{ }

LcfitOnlineAddSequenceMove::~LcfitOnlineAddSequenceMove()
{
    const double lcfit_failure_rate = static_cast<double>(lcfit_failures_) / lcfit_attempts_;
    std::clog << "[LcfitOnlineAddSequenceMove] lcfit failure rate = "
              << lcfit_failures_ << "/" << lcfit_attempts_
              << " (" << lcfit_failure_rate * 100.0 << "%)\n";
}

const std::tuple<bpp::Node*, double, double, double>
LcfitOnlineAddSequenceMove::chooseMoveLocation(bpp::TreeTemplate<bpp::Node>& tree,
                                               const std::string& leafName,
                                               smc::rng* rng,
                                               size_t particleID)
{
  bpp::Node* n = nullptr;
  double edgeLogDensity;
  //std::tie(n, edgeLogDensity) = chooseEdge(tree, leafName, rng, particleID);
    size_t toAddCount = std::distance(taxaToAdd.begin(),taxaToAdd.end());
    
    if( _toAddCount == toAddCount && _probs.find(particleID) != _probs.end() ){
        std::vector<bpp::Node*> nodes = onlineAvailableEdges(tree);
        
        const std::vector<std::pair<size_t, double>>& probabilities = _probs[particleID];
        WeightedSelector<bpp::Node*> selector{*rng};
        for(bpp::Node* node : nodes){
            auto it = std::find_if( probabilities.begin(), probabilities.end(),
                                   [node](const std::pair<size_t, double>& element){return element.first == node->getId();});
            selector.push_back(node, it->second);
        }
        n = selector.choice();
        auto it = std::find_if( probabilities.begin(), probabilities.end(),
                               [n](const std::pair<size_t, double>& element){return element.first == n->getId();});
        edgeLogDensity = log(it->second);
    }
    else{
        if(_toAddCount != toAddCount){
            _probs.clear();
        }
        std::tie(n, edgeLogDensity) = chooseEdge(tree, leafName, rng, particleID);
    }
    
   _toAddCount = toAddCount;
   assert(n);
   assert(std::isfinite(edgeLogDensity));
    
   double mlPendant = 0;
   double mlDistal = 0;
   double distalBranchLength;
   double distalLogDensity;
    
   optimizeBranchLengths(n, leafName, mlDistal, mlPendant);
   std::tie(distalBranchLength, distalLogDensity) = proposeDistal(n->getDistanceToFather(), mlDistal, rng);

   return std::make_tuple(n, edgeLogDensity, distalBranchLength, distalLogDensity);
}

double attachment_lnl_callback(double t, void* data)
{
    AttachmentLikelihood* al = static_cast<AttachmentLikelihood*>(data);

    return (*al)(t);
}

AttachmentProposal LcfitOnlineAddSequenceMove::propose(const std::string& leafName, smc::particle<TreeParticle>& particle, smc::rng* rng)
{
    TreeParticle* treeParticle = particle.GetValuePointer();

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

    bpp::Node* n = nullptr;
    double edgeLogDensity = 0.0;
    double distalBranchLength = 0.0;
    double distalLogDensity = 0.0;

    std::tie(n, edgeLogDensity, distalBranchLength, distalLogDensity) = chooseMoveLocation(*(treeParticle->tree), leafName, rng, treeParticle->particleID);

    assert(n);
    assert(std::isfinite(distalBranchLength));
    assert(std::isfinite(distalLogDensity));

    AttachmentLikelihood al(calculator, n, leafName, distalBranchLength);

    // FIXME: Are there actual branch length constraints available somewhere?
    const double min_t = 1e-6;
    const double max_t = 20.0;
    bsm_t model = DEFAULT_INIT;

    lcfit_fit_auto(&attachment_lnl_callback, &al, &model, min_t, max_t);

    const double mlPendantBranchLength = lcfit_bsm_ml_t(&model);
    double pendantBranchLength, pendantLogDensity;

    try {
        ++lcfit_attempts_;
        lcfit::rejection_sampler sampler(rng->GetRaw(), model, 1.0 / expPriorMean_);
        pendantBranchLength = sampler.sample();
        pendantLogDensity = sampler.log_density(pendantBranchLength);
    } catch (const std::exception& e) {
        // std::clog << "** " << e.what() << '\n';
        ++lcfit_failures_;

        // Fall back on original proposal
        return GuidedOnlineAddSequenceMove::propose(leafName, particle, rng);
    }

    assert(std::isfinite(pendantBranchLength));
    assert(std::isfinite(pendantLogDensity));

    return AttachmentProposal { n, edgeLogDensity, distalBranchLength, distalLogDensity, pendantBranchLength, pendantLogDensity, -1.0, mlPendantBranchLength, "LcfitOnlineAddSequenceMove" };
}

}} // namespace sts::online

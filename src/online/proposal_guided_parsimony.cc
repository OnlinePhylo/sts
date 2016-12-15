#include "proposal_guided_parsimony.h"

#include <limits>


#include "tripod_optimizer.h"
#include "online_util.h"
#include "tree_particle.h"
#include "weighted_selector.h"
#include "composite_tree_likelihood.h"
#include "attachment_likelihood.h"

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

#include <lcfit_rejection_sampler.h>
#include <lcfit_select.h>

namespace sts {
    namespace online {
        
        ProposalGuidedParsimony::~ProposalGuidedParsimony()
        {
            const double lcfit_failure_rate = static_cast<double>(lcfit_failures_) / lcfit_attempts_;
            std::clog << "[ProposalGuidedParsimony] lcfit failure rate = "
            << lcfit_failures_ << "/" << lcfit_attempts_
            << " (" << lcfit_failure_rate * 100.0 << "%)\n";
        }
        
        double attachment_lnl_callback2(double t, void* data)
        {
            AttachmentLikelihood* al = static_cast<AttachmentLikelihood*>(data);
            
            return (*al)(t);
        }
        
        AttachmentProposal ProposalGuidedParsimony::propose(const std::string& leafName, smc::particle<TreeParticle>& particle, smc::rng* rng){
            TreeParticle* value = particle.GetValuePointer();
            std::unique_ptr<bpp::TreeTemplate<bpp::Node>>& tree = value->tree;
            
            bpp::Node* n = nullptr;
            double edgeLogDensity;
            
            size_t toAddCount = std::distance(taxaToAdd.begin(),taxaToAdd.end());
            
            if(_toAddCount == toAddCount && _probs.find(value->particleID) != _probs.end() ){
                std::vector<bpp::Node*> nodes = onlineAvailableEdges(*tree);
                
                std::vector<std::string> names = _parsimony->getNames();
                size_t counter = names.size();
                for(bpp::Node* node : tree->getNodes()){
                    if(node->isLeaf()){
                        size_t pos = find(names.begin(), names.end(), node->getName()) - names.begin();
                        node->setId(static_cast<int>(pos));
                    }
                    else {
                        node->setId(static_cast<int>(counter));
                        counter++;
                    }
                }
                 
                const std::vector<std::pair<size_t, double>>& probabilities = _probs[value->particleID];
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
                std::tie(n, edgeLogDensity) = chooseEdge(*tree, leafName, rng, value->particleID);
            }
            
            double pendantBranchLength;
            double pendantLogDensity;
            double distal;
            double distalLogDensity = 0.0;
            
            double mlDistal = 0;
            double mlPendant = 0;
            
            calculator.calculateAttachmentLikelihood(leafName, n, 0, {0.0});
            
            optimizeBranchLengths(n, leafName, mlDistal, mlPendant);
            
            std::tie(distal, distalLogDensity) = proposeDistal(n->getDistanceToFather(), mlDistal, rng);
            
            AttachmentLikelihood al(calculator, n, leafName, distal);
            
            // FIXME: Are there actual branch length constraints available somewhere?
            const double min_t = 1e-6;
            const double max_t = 20.0;
            bsm_t model = DEFAULT_INIT;
            
            lcfit_fit_auto(&attachment_lnl_callback2, &al, &model, min_t, max_t);
            
            try {
                ++lcfit_attempts_;
                lcfit::rejection_sampler sampler(rng->GetRaw(), model, 1.0 / 0.1);
                pendantBranchLength = sampler.sample();
                pendantLogDensity = sampler.log_density(pendantBranchLength);
            } catch (const std::exception& e) {
                // std::clog << "** " << e.what() << '\n';
                ++lcfit_failures_;
                
                // Fall back on original proposal
                std::tie(pendantBranchLength, pendantLogDensity) = proposePendant(mlPendant, rng);
            }
            
            
            return AttachmentProposal { n, edgeLogDensity, distal, distalLogDensity, pendantBranchLength, pendantLogDensity, mlDistal, mlPendant, "ParsimonyGuidedOnlineAddSequenceMove" };
        }
        
        const std::pair<bpp::Node*, double> ProposalGuidedParsimony::chooseEdge(bpp::TreeTemplate<bpp::Node>& tree, const std::string& leafName, smc::rng* rng, size_t particleID) {
            
            std::vector<bpp::Node*> nodes = onlineAvailableEdges(tree);
            
            if(nodes.size() == 1 ){
                return std::pair<bpp::Node*,double>(nodes[0], 0);
            }
            
            std::vector<bpp::Node*> r = tree.getNodes();
            std::vector<std::string> names = _parsimony->getNames();
            size_t counter = names.size();
            for(bpp::Node* node : r){
                if(node->isLeaf()){
                    size_t pos = find(names.begin(), names.end(), node->getName()) - names.begin();
                    node->setId(static_cast<int>(pos));
                }
                else {
                    node->setId(static_cast<int>(counter));
                    counter++;
                }
            }
            
            std::vector<std::pair<bpp::Node*, double> > nodeWeights;
            double minWeight = std::numeric_limits<double>::max();
            
            double score = _parsimony->getScore(tree);
            
            for(bpp::Node *node : nodes){
                double score = _parsimony->getScore(tree, *node, leafName);
                nodeWeights.push_back(std::make_pair(node, score));
                if(score < minWeight){
                    minWeight = score;
                }
            }
            std::vector<std::pair<bpp::Node*, double> > vec;
            
            double scaler = 0.05;
            double sum = 0.0;
            for(auto &v : nodeWeights){
                double p = exp(scaler*(minWeight-v.second));
                sum += p;
                vec.push_back(std::make_pair(v.first, p));
            }
            
            std::vector<std::pair<size_t, double>> probabilities;
            probabilities.reserve(vec.size());
            
            WeightedSelector<bpp::Node*> selector{*rng};
            for(int i = 0; i < vec.size(); i++){
                double p = vec[i].second/sum;
                selector.push_back(vec[i].first, p);
                probabilities.push_back(std::make_pair(vec[i].first->getId(), p));
            }
            
            _probs[particleID] = probabilities;
            
            bpp::Node* n = selector.choice();
            auto it = std::find_if( nodeWeights.begin(), nodeWeights.end(),
                                   [n](const std::pair<bpp::Node*, double>& element){return element.first == n;});
            double l = scaler*(minWeight-it->second) - log(sum);
            
            return std::pair<bpp::Node*,double>(n, l);
        }
    }
}

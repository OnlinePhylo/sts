#include "proposal_guided_parsimony.h"

#include <limits>


#include "tripod_optimizer.h"
#include "online_util.h"
#include "tree_particle.h"
#include "weighted_selector.h"
#include "composite_tree_likelihood.h"

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

namespace sts {
    namespace online {
        
        /// Propose branch-lengths around ML value
        TripodOptimizer ProposalGuidedParsimony::optimizeBranchLengths(const bpp::Node* insertEdge,
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
            
            return optim;
        }
        
        /// Distal branch length proposal
        /// We propose from Gaussian(mlDistal, edgeLength / 4) with support truncated to [0, edgeLength]
        std::pair<double, double> ProposalGuidedParsimony::proposeDistal(const double edgeLength, const double mlDistal, smc::rng* rng) const
        {
            assert(mlDistal <= edgeLength);
            // HACK/ARBITRARY: proposal standard deviation
            const double sigma = edgeLength / 4;
            
            double distal = -1;
            
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
        
        std::pair<double, double> ProposalGuidedParsimony::proposePendant(const double mlPendant, smc::rng* rng) const
        {
            const double pendantBranchLength = rng->Exponential(mlPendant);
            const double pendantLogDensity = std::log(gsl_ran_exponential_pdf(pendantBranchLength, mlPendant));
            return std::pair<double, double>(pendantBranchLength, pendantLogDensity);
        }
        
        AttachmentProposal ProposalGuidedParsimony::propose(const std::string& leafName, smc::particle<TreeParticle>& particle, smc::rng* rng){
            TreeParticle* value = particle.GetValuePointer();
            std::unique_ptr<bpp::TreeTemplate<bpp::Node>>& tree = value->tree;
            
            bpp::Node* n = nullptr;
            double edgeLogDensity;
            std::tie(n, edgeLogDensity) = chooseEdge(*tree, leafName, rng);
            
            double pendantBranchLength;
            double pendantLogDensity;
            double distal;
            double distalLogDensity = 0.0;
            
            double mlDistal = 0;
            double mlPendant = 0;
            
            if(!_hybrid){
                std::tie(pendantBranchLength, pendantLogDensity) = _branchLengthProposer(rng);
                distal = rng->UniformS() * n->getDistanceToFather();
            }
            else{
                optimizeBranchLengths(n, leafName, mlDistal, mlPendant);
                
                std::tie(distal, distalLogDensity) = proposeDistal(n->getDistanceToFather(), mlDistal, rng);
                
                std::tie(pendantBranchLength, pendantLogDensity) = proposePendant(mlPendant, rng);
            }
            return AttachmentProposal { n, edgeLogDensity, distal, distalLogDensity, pendantBranchLength, pendantLogDensity, mlDistal, mlPendant, "ParsimonyGuidedOnlineAddSequenceMove" };
        }
        
        const std::pair<bpp::Node*, double> ProposalGuidedParsimony::chooseEdge(bpp::TreeTemplate<bpp::Node>& tree, const std::string& leafName, smc::rng* rng) const{
            
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
            
            WeightedSelector<bpp::Node*> selector{*rng};
            for(int i = 0; i < vec.size(); i++){
                selector.push_back(vec[i].first, vec[i].second/sum);
            }
            bpp::Node* n = selector.choice();
            auto it = std::find_if( nodeWeights.begin(), nodeWeights.end(),
                                   [n](const std::pair<bpp::Node*, double>& element){return element.first == n;});
            double l = scaler*(minWeight-it->second) - log(sum);
            
            return std::pair<bpp::Node*,double>(n, l);
        }
    }
}

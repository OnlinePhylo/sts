#include "proposal_guided_parsimony.h"

#include <limits>


#include "online_util.h"
#include "weighted_selector.h"

namespace sts {
    namespace online {
        
        const std::pair<bpp::Node*, double> ProposalGuidedParsimony::chooseEdge(bpp::TreeTemplate<bpp::Node>& tree, const std::string& leafName, smc::rng* rng, size_t particleID) {
            
            std::vector<bpp::Node*> nodes = onlineAvailableEdges(tree);
            
            if(nodes.size() == 1 ){
                return std::pair<bpp::Node*,double>(nodes[0], 0);
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

            double sum = 0.0;
            for(auto &v : nodeWeights){
                double p = exp(_heating*(minWeight-v.second));
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
            double l = _heating*(minWeight-it->second) - log(sum);
            
            return std::pair<bpp::Node*,double>(n, l);
        }
    }
}

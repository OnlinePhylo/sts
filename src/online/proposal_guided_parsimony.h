#ifndef PROPOSAL_GUIDED_PARSIMONY_H
#define PROPOSAL_GUIDED_PARSIMONY_H

#include <stdio.h>

#include "online_add_sequence_move.h"
#include "tripod_optimizer.h"
#include "flexible_parsimony.h"

namespace sts {
    namespace online {
        class ProposalGuidedParsimony : public OnlineAddSequenceMove{
        
        public:
            ProposalGuidedParsimony(std::shared_ptr<FlexibleParsimony> parsimony,
                                    CompositeTreeLikelihood& calculator,
                                    const std::vector<std::string>& taxaToAdd,
                                    std::function<std::pair<double,double>(smc::rng*)> branchLengthProposer, bool hybrid=true): OnlineAddSequenceMove(calculator, taxaToAdd), _parsimony(parsimony), _branchLengthProposer(branchLengthProposer), _hybrid(hybrid){}
            
            virtual AttachmentProposal propose(const std::string& leafName, smc::particle<TreeParticle>& particle, smc::rng* rng);
            
        protected:
            virtual TripodOptimizer optimizeBranchLengths(const bpp::Node* insertEdge,
                                                                           const std::string& newLeafName,
                                                                           double& distalBranchLength,
                                                                           double& pendantBranchLength);
            
            std::pair<double, double> proposeDistal(const double edgeLength, const double mlDistal, smc::rng* rng) const;
            
            std::pair<double, double> proposePendant(const double mlPendant, smc::rng* rng) const;
            
            
        private:
            
            const std::pair<bpp::Node*, double> chooseEdge(bpp::TreeTemplate<bpp::Node>& tree, const std::string& leafName, smc::rng* rng) const;
            
            std::shared_ptr<FlexibleParsimony> _parsimony;
            
            std::function<std::pair<double,double>(smc::rng*)> _branchLengthProposer;
            
            bool _hybrid;
        };
    }
}
#endif

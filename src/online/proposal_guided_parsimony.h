#ifndef PROPOSAL_GUIDED_PARSIMONY_H
#define PROPOSAL_GUIDED_PARSIMONY_H

#include <stdio.h>

#include "guided_online_add_sequence_move.h"
#include "flexible_parsimony.h"

namespace sts {
    namespace online {
        class ProposalGuidedParsimony : public GuidedOnlineAddSequenceMove{
        
        public:
            ProposalGuidedParsimony(std::shared_ptr<FlexibleParsimony> parsimony,
                                    CompositeTreeLikelihood& calculator,
                                    const std::vector<std::string>& taxaToAdd):
            GuidedOnlineAddSequenceMove(calculator, taxaToAdd), _parsimony(parsimony), lcfit_failures_(0), lcfit_attempts_(0){
                _toAddCount = -1;}
            
            virtual ~ProposalGuidedParsimony();
            
            virtual AttachmentProposal propose(const std::string& leafName, smc::particle<TreeParticle>& particle, smc::rng* rng);
            
            virtual const std::pair<bpp::Node*, double> chooseEdge(bpp::TreeTemplate<bpp::Node>& tree,
                                                                   const std::string& leafName,
                                                                   smc::rng* rng, size_t particleID);
            
        private:
            
            std::shared_ptr<FlexibleParsimony> _parsimony;
            
            size_t lcfit_failures_;
            size_t lcfit_attempts_;
        };
    }
}
#endif

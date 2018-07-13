#ifndef PROPOSAL_GUIDED_PARSIMONY_H
#define PROPOSAL_GUIDED_PARSIMONY_H

#include <stdio.h>

#include "lcfit_online_add_sequence_move.h"
#include "flexible_parsimony.h"

namespace sts {
    namespace online {
        class ProposalGuidedParsimony : public LcfitOnlineAddSequenceMove{
        
        public:
            ProposalGuidedParsimony(std::shared_ptr<FlexibleParsimony> parsimony,
                                    std::vector<std::unique_ptr<CompositeTreeLikelihood>>& calculator,
                                    const std::vector<std::string>& sequenceNames,
                                    const std::vector<std::string>& taxaToAdd,
                                    const double expPriorMean):
            LcfitOnlineAddSequenceMove(calculator, sequenceNames, taxaToAdd, {0}, std::numeric_limits<double>::max(),0,expPriorMean), _parsimony(parsimony){
                _proposalMethodName = "ParsimonyOnlineAddSequenceMove";}
            
            virtual ~ProposalGuidedParsimony(){};
            
            virtual const std::pair<bpp::Node*, double> chooseEdge(bpp::TreeTemplate<bpp::Node>& tree,
                                                                   const std::string& leafName,
                                                                   smc::rng* rng, size_t particleID);
            
        private:
            
            std::shared_ptr<FlexibleParsimony> _parsimony;
        };
    }
}
#endif

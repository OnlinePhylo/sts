#ifndef FlexibleTreeLikelihood_hpp
#define FlexibleTreeLikelihood_hpp

#include <stdio.h>
#include <vector>

#include "likelihood.h"

#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/SitePatterns.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>

namespace sts {
    namespace online {
        
        class FlexibleTreeLikelihood: public Likelihood{
            
        public:
            
            virtual ~FlexibleTreeLikelihood(){}
            
            virtual void initialize(const bpp::SubstitutionModel &model, const bpp::DiscreteDistribution& rateDist, bpp::TreeTemplate<bpp::Node>& tree) = 0;
            
            virtual double calculateLogLikelihood() = 0;
            
            virtual double calculateLogLikelihood(const bpp::Node& distal, std::string taxonName, double pendantLength, double distalLength, double proximalLength) = 0;
            
            // Compute derivatives at node node with current tree
//            virtual void calculateDerivatives(const bpp::Node& node, double* d1, double* d2) = 0;
            
            // Compute derivatives of pendant branch with taxon taxonName
            virtual void calculatePendantDerivatives(const bpp::Node& distal, std::string taxonName, double pendantLength, double distalLength, double proximalLength, double* d1, double* d2) = 0;
            
            virtual void calculateDistalDerivatives(const bpp::Node& distal, std::string taxonName, double pendantLength, double distalLength, double proximalLength, double* d1, double* d2) = 0;
			
			virtual void updateNode(const bpp::Node& node) = 0;
			
			virtual void updateAllNodes() = 0;
        };
    }
}

#endif /* FlexibleTreeLikelihood_hpp */

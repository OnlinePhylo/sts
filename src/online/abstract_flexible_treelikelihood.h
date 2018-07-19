#ifndef AbstractFlexibleTreeLikelihood_hpp
#define AbstractFlexibleTreeLikelihood_hpp

#include <stdio.h>
#include <vector>

#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/SitePatterns.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>

#include "flexible_tree_likelihood.h"

namespace sts {
    namespace online {
        
        class AbstractFlexibleTreeLikelihood : public FlexibleTreeLikelihood{
            
        public:
            AbstractFlexibleTreeLikelihood(const bpp::SitePatterns& patterns, const bpp::SubstitutionModel &model, const bpp::DiscreteDistribution& rateDist, bool useAmbiguities=false);
            
            virtual ~AbstractFlexibleTreeLikelihood(){}
            
            virtual void initialize(const bpp::SubstitutionModel &model,const  bpp::DiscreteDistribution& rateDist, bpp::TreeTemplate<bpp::Node>& tree);
            
            virtual double calculateLogLikelihood() = 0;
            
            virtual double calculateLogLikelihood(const bpp::Node& distal, std::string taxonName, double pendantLength, double distalLength, double proximalLength) = 0;
            
            // Compute derivatives at node node with current tree
//            virtual void calculateDerivatives(const bpp::Node& node, double* d1, double* d2) = 0;
            
            // Compute derivatives of pendant branch with taxon taxonName
            virtual void calculatePendantDerivatives(const bpp::Node& distal, std::string taxonName, double pendantLength, double distalLength, double proximalLength, double* d1, double* d2) = 0;
            
            virtual void calculateDistalDerivatives(const bpp::Node& distal, std::string taxonName, double pendantLength, double distalLength, double proximalLength, double* d1, double* d2) = 0;
			
			virtual void updateNode(const bpp::Node& node);
			
			virtual void updateAllNodes();
			
        public:
            static size_t operationCallCount;
            
        protected:
	        
	        virtual void setStates(const bpp::SiteContainer& sites) = 0;
	        
	        virtual void setPartials(const bpp::SiteContainer& sites) = 0;
	        
	        
            bpp::SubstitutionModel const* _model;
            bpp::DiscreteDistribution const* _rateDist;
            bpp::TreeTemplate<bpp::Node>* _tree;
            
            std::vector<std::string> _taxa;
            
            bool _useAmbiguities;
            
            int _stateCount;
            int _patternCount;
            int _sequenceCount;
            int _totalNodeCount; // number of nodes in full tree
            int _rateCount;
            
            int _partialCount; // includes upperPartials
            int _matrixCount; // temporary matrices
            
            size_t _nodeCount; // current number of nodes
            size_t _leafCount; // current number of leaves
            size_t _internalNodeCount;
            
            std::vector<double> _patternWeights;
            
            bool _updateSubstitutionModel;
            bool _updateSiteModel;
            std::vector<bool> _needNodeUpdate;
            bool _updatePartials;
            bool _updateUpperPartials;
            
            bool _useScaleFactors;
            bool _recomputeScaleFactors;
        };
    }
}
    
#endif /* AbstractFlexibleTreeLikelihood_hpp */

#ifndef BeagleFlexibleTreeLikleihood_hpp
#define BeagleFlexibleTreeLikleihood_hpp

#include <stdio.h>
#include <vector>


#include "abstract_flexible_treelikelihood.h"

#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/SitePatterns.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>

#include "libhmsbeagle/beagle.h"

namespace sts {
    namespace online {
        
        class BeagleFlexibleTreeLikelihood : public AbstractFlexibleTreeLikelihood {
            
        public:
            BeagleFlexibleTreeLikelihood(const bpp::SitePatterns& patterns, const bpp::SubstitutionModel &model, const bpp::DiscreteDistribution& rateDist, bool useAmbiguities=true);
            
            virtual ~BeagleFlexibleTreeLikelihood();
            
            virtual void initialize(const bpp::SubstitutionModel &model, const bpp::DiscreteDistribution& rateDist, bpp::TreeTemplate<bpp::Node>& tree);
            
            virtual double calculateLogLikelihood();
            
            virtual double calculateLogLikelihood(double length);
            
            virtual double calculateLogLikelihood(const bpp::Node& distal, std::string taxonName, double pendantLength, double distalLength, double proximalLength);
            
//            virtual void calculateDerivatives(const bpp::Node& node, double* d1, double* d2);
            
            // Compute derivatives of pendant branch with taxon taxonName
            virtual void calculatePendantDerivatives(const bpp::Node& distal, std::string taxonName, double pendantLength, double distalLength, double proximalLength, double* d1, double* d2);
            
            virtual void calculateDistalDerivatives(const bpp::Node& distal, std::string taxonName, double pendantLength, double distalLength, double proximalLength, double* d1, double* d2);
            
        protected:
            
            void updateSubstitutionModel();
            
            void updateSiteModel();
            
            bool traverse(const bpp::Node* node);
            
            void traverseUpper(const bpp::Node* node);
            
            void setStates(const bpp::SiteContainer& sites);
            
            void setPartials(const bpp::SiteContainer& sites);
            
        private:
            double _logLnl;
            
            
            int _beagleInstance;
            BeagleInstanceDetails _instanceDetails;
            
            //size_t _nBeagleUpdateTransitionsCalls;
            
            std::vector<int> _matrixUpdateIndices;
            std::vector<double> _branchLengths;
            std::vector<int> _scaleBufferIndices;
            std::vector<int> _scaleBufferUpperIndices;
            std::vector<BeagleOperation> _operations;

            std::vector<int> _upperPartialsIndexes;
            
        };
    }
}

#endif /* BeagleFlexibleTreeLikleihood_hpp */

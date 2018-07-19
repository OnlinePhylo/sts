#ifndef SimpleFlexibleTreeLikelihood_hpp
#define SimpleFlexibleTreeLikelihood_hpp

#include <stdio.h>
#include <vector>

#include "abstract_flexible_treelikelihood.h"

#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/SitePatterns.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>

namespace sts {
    namespace online {
    
        class SimpleFlexibleTreeLikelihood : public AbstractFlexibleTreeLikelihood{
            
        public:
            SimpleFlexibleTreeLikelihood(const bpp::SitePatterns& patterns, const bpp::SubstitutionModel &model, const bpp::DiscreteDistribution& rateDist, bool useAmbiguities=true);
            
            virtual ~SimpleFlexibleTreeLikelihood(){}
            
            virtual double calculateLogLikelihood();
            
            virtual double calculateLogLikelihood(const bpp::Node& distal, std::string taxonName, double pendantLength, double distalLength, double proximalLength);
            
            virtual void calculatePendantDerivatives(const bpp::Node& distal, std::string taxonName, double pendantLength, double distalLength, double proximalLength, double* d1, double* d2);

            virtual void calculateDistalDerivatives(const bpp::Node& distal, std::string taxonName, double pendantLength, double distalLength, double proximalLength, double* d1, double* d2);
    
        protected:
            
            virtual void setStates(const bpp::SiteContainer& sites);
            
            virtual void setPartials(const bpp::SiteContainer& sites);
            
            bool traverse(const bpp::Node* node);
            
            void traverseUpper(const bpp::Node* node);
            
            void calculateBranchLikelihood(double* rootPartials, const double* attachmentPartials, const double* pendantPartials, const double* pendantMatrices, const double* weights);
            
            void calculateBranchLikelihood(double* rootPartials, const double* attachmentPartials, const int* pendantStates, const double* pendantMatrices, const double* weights);
                
            
            
            void calculatePatternLikelihood( const double *partials, const double *frequencies, double *outLogLikelihoods)const;
            
            void integratePartials( const double *inPartials, const double *proportions, double *outPartials )const;
            
            void updatePartialsKnownKnown( const int *states1, const double *matrices1, const int *states2, const double *matrices2, double *partials )const;
            
            void updatePartialsKnownUndefined( const int *states1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 )const;
            
            void updatePartialsUndefinedUndefined( const double *partials1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 )const;            
            
            void updatePartials(int partialsIndex, int partialsIndex1, int matrixIndex1, int partialsIndex2, int matrixIndex2 );
            
            
//            void updateUpperPartialsKnown( const double *matrix_upper, const double *partials_upper, const double *matrix_lower, const int *states, double *partials ) const;
//            
//            void updateUpperPartialsUndefined( const double *matrix_upper, const double *partials_upper, const double *matrix_lower, const double *partials_lower, double *partials ) const;
//            
//            void updateUpperPartialsRootUndefined( const double *partials1, const double *matrices1, const double *frequencies, double *partials_upper )const;
//            
//            void updateUpperPartialsRootKnown(const double *matrices1, const int *states1, const double *frequencies, double *partials_upper ) const;
//            
//            void updateUpperPartials(const bpp::Node* node);
//            
//            void calculatePatternLikelihood( const double *partials_upper, const double *partials_lower, const double *matrix_lower, const double *proportions, double *pattern_lk ) const;
//            
//            void calculatePatternLikelihoodLeaf( const double *partials_upper, const int* states, const double *matrix_lower, const double *proportions, double *pattern_lk ) const;
            
        private:
//            void calculateDerivatives(const bpp::Node& distal, std::string taxonName, int index, double pendantLength, double distalLength, double proximalLength, double* d1, double* d2);
            
            std::vector<std::vector<int> > _states;

            double _logLnl;
            
            int _matrixSize;
            
            std::vector<double> _patternWeights;
            
            std::vector<std::vector<double> > _matrices;
            std::vector<std::vector<double> > _partials;
            std::vector<std::vector<double> > _rootPartials;
            std::vector<std::vector<double> > _patternLikelihoods;
            
            std::vector<int> _upperPartialsIndexes;
        };
    }
}
#endif /* SimpleFlexibleTreeLikelihood_hpp */

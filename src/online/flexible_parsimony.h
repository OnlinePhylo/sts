#ifndef FLEXIBLE_PARSIMONY_H
#define FLEXIBLE_PARSIMONY_H

#include <stdio.h>

#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/SitePatterns.h>

namespace sts {
    namespace online{
        
        class FlexibleParsimony{

        public:
            FlexibleParsimony(const bpp::SitePatterns& patterns, const bpp::Alphabet& alphabet);
            
            FlexibleParsimony() = delete;
            
            FlexibleParsimony(const FlexibleParsimony&) = delete;            //disable copy-constructor
            
            FlexibleParsimony& operator=(const FlexibleParsimony&) = delete; //disable copy-assignment
            
            virtual ~FlexibleParsimony(){}
            
            double getScore(const bpp::TreeTemplate<bpp::Node>& tree);
            
            double getScore(const bpp::TreeTemplate<bpp::Node>& tree, const bpp::Node& distal, std::string taxon);
            
            void updateAllNodes();
            
            void updateNode(const bpp::Node& node);
            
        protected:
            
            bool first_pass( const bpp::Node& node );
            
            bool first_pass( const bpp::Node& node, const bpp::Node& distal, size_t indexAttachment, size_t indexTaxon );
            
            void calculateLocalScore(int8_t* states, int32_t* local_scores,
                                                        const int8_t* states1, const int8_t* states2,
                                                        const int32_t* local_scores1, const int32_t* local_scores2,
                                                        bool isLeaf1, bool isLeaf2);
//#ifdef __SSE2__
//            void calculateLocalScoreSSE(int8_t* states, int32_t* local_scores,
//                                                           const int8_t* states1, const int8_t* states2,
//                                                           const int32_t* local_scores1, const int32_t* local_scores2,
//                                                           bool isLeaf1, bool isLeaf2);
//#endif
        private:
            
            size_t _stateCount;
            size_t _patternCount;
            size_t _sequenceCount;
            size_t _nodeCount;
            
            std::vector<std::string> _taxa;
            
            double _score;
            bool _updateScores;
            std::vector<bool> _updateNode;
            
            std::vector<std::vector<int8_t> > _stateSets;
            
            std::vector<int32_t> _weights;
            std::vector<std::vector<int32_t> > _local_scores;

        };
    }
}
#endif

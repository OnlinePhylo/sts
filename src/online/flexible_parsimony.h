#ifndef FLEXIBLE_PARSIMONY_H
#define FLEXIBLE_PARSIMONY_H

#include <stdio.h>

#include <Bpp/Phyl/TreeTemplate.h>

namespace sts {
    namespace online{
        
        class FlexibleParsimony{

        public:
            FlexibleParsimony(const bpp::SiteContainer& sites);
            
            FlexibleParsimony() = delete;
            
            FlexibleParsimony(const FlexibleParsimony&) = delete;            //disable copy-constructor
            
            FlexibleParsimony& operator=(const FlexibleParsimony&) = delete; //disable copy-assignment
            
            virtual ~FlexibleParsimony(){}
            
            double getScore(const bpp::TreeTemplate<bpp::Node>& tree);
            
            double getScore(const bpp::TreeTemplate<bpp::Node>& tree, const bpp::Node& distal, std::string taxon);

            const std::vector<std::string> getNames() const;
            
            void updateAllNodes();
            
            void updateNode(const bpp::Node& node);
            
        protected:
            
            bool first_pass( const bpp::Node& node );
            
            bool first_pass( const bpp::Node& node, const bpp::Node& distal, size_t indexAttachment, size_t indexTaxon );
            
            void calculateLocalScore(std::vector<int8_t>& states, std::vector<int32_t>& local_scores,
                                                        const std::vector<int8_t>& states1, const std::vector<int8_t>& states2,
                                                        const std::vector<int32_t>& local_scores1, const std::vector<int32_t>& local_scores2);
            
            void calculateLocalScore(int8_t* states, int32_t* local_scores,
                                                        const int8_t* states1, const int8_t* states2,
                                                        const int32_t* local_scores1, const int32_t* local_scores2,
                                                        bool isLeaf1, bool isLeaf2);
            
        private:
            
            const bpp::SiteContainer& _sites;
            
            
            size_t _stateCount;
            size_t _patternCount;
            size_t _sequenceCount;
            size_t _nodeCount;
            
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

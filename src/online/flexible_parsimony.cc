#include "flexible_parsimony.h"
#include <cstring>

//#ifdef __SSE2__
//#include <xmmintrin.h> // SSE
//#include <emmintrin.h> // SSE2
//#endif

using namespace std;

namespace sts {
    namespace online {
        
        
        FlexibleParsimony::FlexibleParsimony(const bpp::SitePatterns& patterns, const bpp::Alphabet& alphabet){
            _patternCount = patterns.getWeights().size();
            _stateCount = alphabet.getSize();
            std::unique_ptr<bpp::SiteContainer> sites(patterns.getSites());
            _sequenceCount = sites->getNumberOfSequences();
            _nodeCount = (_sequenceCount * 2) - 1; // number of nodes
            
            const size_t stateSetsCount = _nodeCount*2+2;
            const size_t localScoresCount = _nodeCount*2+2;
            
            _stateSets.resize(stateSetsCount);
            for(auto it = _stateSets.begin(); it != _stateSets.end(); ++it){
                it->resize(_stateCount*_patternCount);
            }
            _local_scores.resize(localScoresCount);// #only internal nodes needs storage
            for (size_t i = _sequenceCount; i < localScoresCount; i++) {
                _local_scores[i].assign(_patternCount, 0);
            }
            
            _updateScores = true;
            _updateUpperScores = true;
            _updateNode.assign(_nodeCount, true);
            
            const std::vector<unsigned int>& weights = patterns.getWeights();
            _weights.reserve(_patternCount);
            auto castit = [](unsigned int w) { return static_cast<size_t>(w); };
            std::transform(weights.begin(), weights.end(), _weights.begin(), castit);

            
            _taxa = sites->getSequencesNames();
            
            // initialize stateSets
            for ( int i = 0; i < _sequenceCount; i++ ) {
                const bpp::Sequence& sequence = sites->getSequence(i);
                
                _stateSets[i].assign(_stateCount*_patternCount, 0);
                for ( int j = 0; j < _patternCount; j++ ) {
                    size_t pattern = alphabet.getStateIndex(sequence.getChar(j))-1;
                    if( pattern < _stateCount ){
                        _stateSets[i][j*_stateCount+pattern ] = 1;
                    }
                    else {
                        for ( int k = 0; k < _stateCount; k++ ) {
                            _stateSets[i][j*_stateCount+k] = 1;
                        }
                    }
                }
                
            }
            _score = 0;
            _upperPartialsIndexes.resize(_nodeCount);
        }
        
        double FlexibleParsimony::getScore(const bpp::TreeTemplate<bpp::Node>& tree){
            _updateScores = true;
            _updateNode.assign(_updateNode.size(), true);
            
            if ( _updateScores ) {
                first_pass(*tree.getRootNode());
                
                _score = 0;
                const std::vector<int32_t>& root_scores =_local_scores[tree.getRootNode()->getId()];
                for ( int i = 0; i < _patternCount; i++ ) {
                    _score += root_scores[i] * _weights[i];
                }
                _updateNode.assign(tree.getNumberOfNodes(), false);
                _updateScores = false;
                _updateUpperScores = true;
            }
            return _score;
        }
        
        double FlexibleParsimony::getScore(const bpp::TreeTemplate<bpp::Node>& tree, const bpp::Node& distal, std::string taxon){
            if(_updateUpperScores){
                traverseUpper(tree.getRootNode());
                _updateUpperScores = false;
            }
            
            size_t indexTaxon = find(_taxa.begin(), _taxa.end(), taxon) - _taxa.begin();
            assert(tree.getRootNode()->getId() != _local_scores.size()-1);
            
            
            // Connect distal and proximal (upper)
            const size_t tempIdx = _local_scores.size()-2;
            const size_t tempRootIdx = _local_scores.size()-1;
            calculateLocalScore(_stateSets[tempIdx].data(), _local_scores[tempIdx].data(),
                                _stateSets[distal.getId()].data(), _stateSets[_upperPartialsIndexes[distal.getId()]].data(),
                                _local_scores[distal.getId()].data(), _local_scores[_upperPartialsIndexes[distal.getId()]].data(),
                                _local_scores[distal.getId()].size() == 0, _local_scores[_upperPartialsIndexes[distal.getId()]].size()==0);
            
            // Connect to pendant
            calculateLocalScore(_stateSets[tempRootIdx].data(), _local_scores[tempRootIdx].data(),
                                _stateSets[tempIdx].data(), _stateSets[indexTaxon].data(),
                                _local_scores[tempIdx].data(), _local_scores[indexTaxon].data(),
                                _local_scores[tempIdx].size() == 0, true);
            
            _score = 0;
            const std::vector<int32_t>& root_scores =_local_scores[tempRootIdx];
            for ( int i = 0; i < _patternCount; i++ ) {
                _score += root_scores[i] * _weights[i];
            }

            return _score;
        }
        
        bool FlexibleParsimony::first_pass( const bpp::Node& node ){
            bool updated = _updateNode[node.getId()];
            
            if( node.getNumberOfSons() > 0 ) {
                bool updated_child = first_pass(*node.getSon(0));
                updated_child |= first_pass(*node.getSon(1));
                
                if( updated_child ){
                    int nodeId = node.getId();
                    size_t idx0 = node.getSon(0)->getId();
                    size_t idx1 = node.getSon(1)->getId();
                    std::vector<int8_t>& states = _stateSets[nodeId];
                    
                    for ( int i = 0; i < _patternCount; i++ ) {
                        _local_scores[nodeId][i] = 0;
                        
                        int size = 0;
                        size_t offset = i*_stateCount;
                        
                        // intersection
                        for ( int j = 0; j < _stateCount; j++ ) {
                            int8_t temp = _stateSets[idx0][offset+j] & _stateSets[idx1][offset+j];
                            states[offset+j] = temp;
                            if( temp ){
                                size++;
                            }
                        }
                        
                        if( size == 0 ){
                            // Union
                            for ( int j = 0; j < _stateCount; j++ ) {
                                states[offset+j] = _stateSets[idx0][offset+j] | _stateSets[idx1][offset+j];
                            }
                            _local_scores[nodeId][i] = 1;
                        }
                        
                    }
                    
                    std::vector<int32_t>& local_score = _local_scores[nodeId];
                    
                    if( !node.getSon(0)->isLeaf() ){
                        for ( int i = 0; i < _patternCount; i++ ){
                            local_score[i] += _local_scores[idx0][i];
                        }
                    }
                    if( !node.getSon(1)->isLeaf() ){
                        for ( int i = 0; i < _patternCount; i++ ){
                            local_score[i] += _local_scores[idx1][i];
                        }
                    }
                    
                    updated = true;
                }
            }
            return updated;
        }

        void FlexibleParsimony::traverseUpper(const bpp::Node* node){
            
            if(node->hasFather()){
                const bpp::Node* parent = node->getFather();
                const bpp::Node* sibling = parent->getSon(0);
                if (sibling == node) {
                    sibling = parent->getSon(1);
                }
                
                if(parent->hasFather()){
                    const bpp::Node* grandParent = parent->getFather();
                    const int idSibling = sibling->getId();
                    
                    int idMatrix = parent->getId();
                    
                    // The sons of the right node of the root are going to use the lower partials of the left node of the root
                    if(!grandParent->hasFather() && grandParent->getSon(1) == parent){
                        idMatrix = grandParent->getSon(0)->getId();
                    }
                    
                    _upperPartialsIndexes[node->getId()] = node->getId() + _nodeCount;
                    
                    calculateLocalScore(_stateSets[_upperPartialsIndexes[node->getId()]].data(), _local_scores[_upperPartialsIndexes[node->getId()]].data(),
                                        _stateSets[_upperPartialsIndexes[parent->getId()]].data(), _stateSets[idSibling].data(),
                                        _local_scores[_upperPartialsIndexes[parent->getId()]].data(), _local_scores[idSibling].data(),
                                        _local_scores[_upperPartialsIndexes[parent->getId()]].size() == 0, _local_scores[idSibling].size() == 0);
                }
                // We dont need to calculate upper partials for the children of the root as it is using the lower partials of its sibling
                // Left node of the root
                else if(parent->getSon(0) == node){
                    _upperPartialsIndexes[node->getId()] = parent->getSon(1)->getId();
                }
                else{
                    _upperPartialsIndexes[node->getId()] = parent->getSon(0)->getId();
                }
            }
            
            if(node->getNumberOfSons() > 0){
                traverseUpper(node->getSon(0));
                traverseUpper(node->getSon(1));
            }
        }
        
        
        void FlexibleParsimony::calculateLocalScore(int8_t* states, int32_t* local_scores,
                                                    const int8_t* states1, const int8_t* states2,
                                                    const int32_t* local_scores1, const int32_t* local_scores2,
                                                    bool isLeaf1, bool isLeaf2){
            
            
            for ( int i = 0; i < _patternCount; i++ ) {
                local_scores[i] = 0;
                
                int size = 0;
                size_t offset = i*_stateCount;
                
                const int8_t* s1 = states1+offset;
                const int8_t* s2 = states2+offset;
                int8_t* s  = states+offset;
                
                // intersection
                for ( int j = 0; j < _stateCount; j++ ) {
                    int8_t temp = *s1 & *s2;
                    *s = temp;
                    if( temp ){
                        size++;
                    }
                    s++; s1++; s2++;
                }
                
                if( size == 0 ){
                    s1 = states1+offset;
                    s2 = states2+offset;
                    s  = states+offset;
                    
                    // Union
                    for ( int j = 0; j < _stateCount; j++ ) {
                        *s = *s1 | *s2;
                        s++; s1++; s2++;
                    }
                    local_scores[i] = 1;
                }
            }
            
            if(!isLeaf1 && !isLeaf2){
                for ( int i = 0; i < _patternCount; i++ ){
                    local_scores[i] += local_scores1[i] + local_scores2[i];
                }
            }
            else if(!isLeaf1){
                for ( int i = 0; i < _patternCount; i++ ) local_scores[i] += local_scores1[i];
            }
            else if(!isLeaf2){
                for ( int i = 0; i < _patternCount; i++ ) local_scores[i] += local_scores2[i];
            }
        }

//#ifdef __SSE2__
//        void FlexibleParsimony::calculateLocalScoreSSE(int8_t* states, int32_t* local_scores,
//                                                    const int8_t* states1, const int8_t* states2,
//                                                    const int32_t* local_scores1, const int32_t* local_scores2,
//                                                    bool isLeaf1, bool isLeaf2){
//            
//            memset(local_scores, 0, sizeof(int32_t)*_patternCount);
//            
//            int i = 0;
//            int n = (_patternCount*4/16)*16;
//            
//            const int8_t* ll = states1;
//            const int8_t* rr = states2;
//            int8_t* s  = states;
//            
//            __m128i *l   = (__m128i*)ll;
//            __m128i *r   = (__m128i*)rr;
//            
//            int8_t temp[16] __attribute__ ((aligned (16)));
//            int32_t temp2[4];
//            
//            int k = 0;
//            
//            for (; i < n; i+=16, l++, r++, s+=16 ) {
//                
//                _mm_store_si128((__m128i*)s, _mm_and_si128(*l, *r));
//                
//                _mm_store_si128((__m128i*)temp, _mm_or_si128(*l, *r));
//                
//                memcpy(temp2, s, sizeof(int8_t)*16);
//                
//                if( temp2[0] == 0 ){
//                    memcpy(s, temp, 4*sizeof(int8_t));
//                    local_scores[k] = 1;
//                }
//                k++;
//                
//                if( temp2[1] == 0 ){
//                    memcpy(s+4, temp+4, 4*sizeof(int8_t));
//                    local_scores[k] = 1;
//                }
//                k++;
//                
//                if( temp2[2] == 0 ){
//                    memcpy(s+8, temp+8, 4*sizeof(int8_t));
//                    local_scores[k] = 1;
//                }
//                k++;
//                
//                if( temp2[3] == 0 ){
//                    memcpy(s+12, temp+12, 4*sizeof(int8_t));
//                    local_scores[k] = 1;
//                }
//                k++;
//            }
//            
//            ll += i;
//            rr += i;
//            
//            for ( i /= 4; i < _patternCount; i++ ) {
//                int size = 0;
//                // intersection
//                s[0] = ll[0] & rr[0];
//                s[1] = ll[1] & rr[1];
//                s[2] = ll[2] & rr[2];
//                s[3] = ll[3] & rr[3];
//                
//                size = s[0] + s[1] + s[2] + s[3];
//                
//                
//                if( size == 0 ){
//                    // Union
//                    s[0] = ll[0] | rr[0];
//                    s[1] = ll[1] | rr[1];
//                    s[2] = ll[2] | rr[2];
//                    s[3] = ll[3] | rr[3];
//                    
//                    local_scores[k] = 1;
//                }
//                s += 4;
//                ll += 4;
//                rr += 4;
//                k++;
//            }
//            
//            __m128i *r1,*r2;
//            
//            if( !isLeaf1){
//                size_t n = (_patternCount/4)*4;
//                r1 = (__m128i *)local_scores;
//                r2 = (__m128i *)local_scores1;
//                size_t i = 0;
//                for ( ; i < n; i+=4 ) {
//                    *r1  = _mm_add_epi32(*r1, *r2);
//                    r1++;r2++;
//                }
//                
//                for ( ; i < _patternCount; i++ ) {
//                    local_scores[i] += local_scores1[i];
//                }
//            }
//            if( !isLeaf2){
//                size_t n = (_patternCount/4)*4;
//                r1 = (__m128i *)local_scores;
//                r2 = (__m128i *)local_scores2;
//                size_t i = 0;
//                for ( ; i < n; i+=4 ) {
//                    *r1  = _mm_add_epi32(*r1, *r2);
//                    r1++;r2++;
//                }
//                
//                for ( ; i < _patternCount; i++ ) {
//                    local_scores[i] += local_scores2[i];
//                }
//            }
//        }
//#endif
        
        bool FlexibleParsimony::first_pass( const bpp::Node& node, const bpp::Node& distal, size_t indexAttachment, size_t indexTaxon ){
            bool updated = _updateNode[node.getId()];
            
            if( node.getNumberOfSons() > 0 ) {
                bool updated_child = first_pass(*node.getSon(0), distal, indexAttachment, indexTaxon);
                updated_child |= first_pass(*node.getSon(1), distal, indexAttachment, indexTaxon);
                
                
                if( updated_child ){
                    int nodeId = node.getId();
                    size_t childCount = node.getNumberOfSons();
                    
                    size_t index1 = node.getSon(0)->getId();
                    size_t index2 = node.getSon(1)->getId();
                    for ( int k = 0; k < childCount; k++ ) {
                        if(*node.getSon(k) == distal){
                            calculateLocalScore(_stateSets[indexAttachment].data(), _local_scores[indexAttachment].data(), _stateSets[distal.getId()].data(), _stateSets[indexTaxon].data(), _local_scores[distal.getId()].data(), _local_scores[indexTaxon].data(), _local_scores[distal.getId()].size() == 0, _local_scores[indexTaxon].size() == 0);
                            index1 = indexAttachment;
                            index2 = node.getSon(1-k)->getId();
                        }
                    }
                    calculateLocalScore(_stateSets[nodeId].data(), _local_scores[nodeId].data(), _stateSets[index1].data(), _stateSets[index2].data(), _local_scores[index1].data(), _local_scores[index2].data(), _local_scores[index1].size() == 0, _local_scores[index2].size() == 0);
                    
                    updated = true;
                }
            }
            return updated;
        }

        void FlexibleParsimony::updateAllNodes(){
            _updateNode.assign(_updateNode.size(), true);
            _updateScores = true;
        }
        
        void FlexibleParsimony::updateNode(const bpp::Node& node){
            _updateNode[node.getId()] = true;
            _updateScores = true;
        }
    }
}

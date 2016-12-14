#include "flexible_parsimony.h"
#include <cstring>

//#ifdef __SSE2__
#include <xmmintrin.h> // SSE
#include <emmintrin.h> // SSE2

//#ifdef __SSE4_1__
#include <smmintrin.h> //SSE4.1
//#endif

using namespace std;

namespace sts {
    namespace online {
        
        
        FlexibleParsimony::FlexibleParsimony(const bpp::SiteContainer& sites):_sites(sites){
            _patternCount = sites.getNumberOfSites();
//            cout<<sites.getAlphabet()->getNumberOfStates()<<endl;
            _stateCount = 4;//sites.getAlphabet()->getNumberOfStates();
            _sequenceCount = sites.getNumberOfSequences();
            _nodeCount = (_sequenceCount * 2) - 1; // number of nodes
            
            _stateSets.resize(_nodeCount);
            for(auto it = _stateSets.begin(); it != _stateSets.end(); ++it){
                it->resize(_stateCount*_patternCount);
            }
            _local_scores.resize(_nodeCount);// #only internal nodes needs storage
            for (size_t i = _sequenceCount; i < _nodeCount; i++) {
                _local_scores[i].assign(_patternCount, 0);
            }
            
            _updateScores = true;
            _updateNode.assign(_nodeCount, true);
            
            _weights.assign(_patternCount, 1.);
            
            // initialize stateSets
            const bpp::Alphabet* alphabet = _sites.getAlphabet();
            for ( int i = 0; i < _sequenceCount; i++ ) {
                const bpp::Sequence& sequence = _sites.getSequence(i);
                
                _stateSets[i].assign(_stateCount*_patternCount, 0);
                for ( int j = 0; j < _patternCount; j++ ) {
                    size_t pattern = alphabet->getStateIndex(sequence.getChar(j))-1;
                    //cout << sequence.getChar(j) << " "<<alphabet->getStateIndex(sequence.getChar(j))<<endl;
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
            }
            return _score;
        }
        
        // Always dirty after this call
        double FlexibleParsimony::getScore(const bpp::TreeTemplate<bpp::Node>& tree, const bpp::Node& distal, std::string taxon){
            //_updateScores = true;
            //_updateNode.assign(_updateNode.size(), true);
            
            _updateScores = true;
            _updateNode[distal.getId()] = true;
            
            size_t indexTaxon = _sites.getSequencePosition(taxon);
            assert(tree.getRootNode()->getId() != _local_scores.size()-1);
            
            first_pass(*tree.getRootNode(), distal, _local_scores.size()-1, indexTaxon);
            
            
            _score = 0;
            const std::vector<int32_t>& root_scores =_local_scores[tree.getRootNode()->getId()];
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

        void FlexibleParsimony::calculateLocalScore(std::vector<int8_t>& states, std::vector<int32_t>& local_scores,
                                                    const std::vector<int8_t>& states1, const std::vector<int8_t>& states2,
                                                    const std::vector<int32_t>& local_scores1, const std::vector<int32_t>& local_scores2){
            
            for ( int i = 0; i < _patternCount; i++ ) {
                local_scores[i] = 0;
                
                int size = 0;
                size_t offset = i*_stateCount;
                
                // intersection
                for ( int j = 0; j < _stateCount; j++ ) {
                    int8_t temp = states1[offset+j] & states2[offset+j];
                    states[offset+j] = temp;
                    if( temp ){
                        size++;
                    }
                }
                
                if( size == 0 ){
                    // Union
                    for ( int j = 0; j < _stateCount; j++ ) {
                        int8_t temp = states1[offset+j] | states2[offset+j];
                        states[offset+j] = temp;
                    }
                    local_scores[i] = 1;
                }
            }
            
            if(local_scores1.size() > 0){
                std::transform(local_scores.begin(), local_scores.end(),
                               local_scores1.begin(),
                               local_scores.begin(),
                               std::plus<int32_t>());
            }
            if(local_scores2.size() > 0){
                std::transform(local_scores.begin(), local_scores.end(),
                               local_scores2.begin(),
                               local_scores.begin(),
                               std::plus<int32_t>());
            }
        }
        
        /*void FlexibleParsimony::calculateLocalScore(int8_t* states, int32_t* local_scores,
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
        }*/
        
        void FlexibleParsimony::calculateLocalScore(int8_t* states, int32_t* local_scores,
                                                    const int8_t* states1, const int8_t* states2,
                                                    const int32_t* local_scores1, const int32_t* local_scores2,
                                                    bool isLeaf1, bool isLeaf2){
            
            memset(local_scores, 0, sizeof(int32_t)*_patternCount);
            
            int i = 0;
            int n = (_patternCount*4/16)*16;
            
            const int8_t* ll = states1;
            const int8_t* rr = states2;
            int8_t* s  = states;
            
            __m128i *l   = (__m128i*)ll;
            __m128i *r   = (__m128i*)rr;
            
            int8_t temp[16] __attribute__ ((aligned (16)));
            int32_t temp2[4];
            
            int k = 0;
            
            for (; i < n; i+=16, l++, r++, s+=16 ) {
                
                _mm_store_si128((__m128i*)s, _mm_and_si128(*l, *r));
                
                _mm_store_si128((__m128i*)temp, _mm_or_si128(*l, *r));
                
                memcpy(temp2, s, sizeof(int8_t)*16);
                
                if( temp2[0] == 0 ){
                    memcpy(s, temp, 4*sizeof(int8_t));
                    local_scores[k] = 1;
                }
                k++;
                
                if( temp2[1] == 0 ){
                    memcpy(s+4, temp+4, 4*sizeof(int8_t));
                    local_scores[k] = 1;
                }
                k++;
                
                if( temp2[2] == 0 ){
                    memcpy(s+8, temp+8, 4*sizeof(int8_t));
                    local_scores[k] = 1;
                }
                k++;
                
                if( temp2[3] == 0 ){
                    memcpy(s+12, temp+12, 4*sizeof(int8_t));
                    local_scores[k] = 1;
                }
                k++;
            }
            
            ll += i;
            rr += i;
            
            for ( i /= 4; i < _patternCount; i++ ) {
                int size = 0;
                // intersection
                s[0] = ll[0] & rr[0];
                s[1] = ll[1] & rr[1];
                s[2] = ll[2] & rr[2];
                s[3] = ll[3] & rr[3];
                
                size = s[0] + s[1] + s[2] + s[3];
                
                
                if( size == 0 ){
                    // Union
                    s[0] = ll[0] | rr[0];
                    s[1] = ll[1] | rr[1];
                    s[2] = ll[2] | rr[2];
                    s[3] = ll[3] | rr[3];
                    
                    local_scores[k] = 1;
                }
                s += 4;
                ll += 4;
                rr += 4;
                k++;
            }
            
            __m128i *r1,*r2;
            
            if( !isLeaf1){
                size_t n = (_patternCount/4)*4;
                r1 = (__m128i *)local_scores;
                r2 = (__m128i *)local_scores1;
                size_t i = 0;
                for ( ; i < n; i+=4 ) {
                    *r1  = _mm_add_epi32(*r1, *r2);
                    r1++;r2++;
                }
                
                for ( ; i < _patternCount; i++ ) {
                    local_scores[i] += local_scores1[i];
                }
            }
            if( !isLeaf2){
                size_t n = (_patternCount/4)*4;
                r1 = (__m128i *)local_scores;
                r2 = (__m128i *)local_scores2;
                size_t i = 0;
                for ( ; i < n; i+=4 ) {
                    *r1  = _mm_add_epi32(*r1, *r2);
                    r1++;r2++;
                }
                
                for ( ; i < _patternCount; i++ ) {
                    local_scores[i] += local_scores2[i];
                }
            }
        }
        
        
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
//                            calculateLocalScore(_stateSets[indexAttachment], _local_scores[indexAttachment], _stateSets[distal.getId()], _stateSets[indexTaxon], _local_scores[distal.getId()], _local_scores[indexTaxon]);
                            calculateLocalScore(_stateSets[indexAttachment].data(), _local_scores[indexAttachment].data(), _stateSets[distal.getId()].data(), _stateSets[indexTaxon].data(), _local_scores[distal.getId()].data(), _local_scores[indexTaxon].data(), _local_scores[distal.getId()].size() == 0, _local_scores[indexTaxon].size() == 0);
                            index1 = indexAttachment;
                            index2 = node.getSon(1-k)->getId();
                        }
                    }
//                    calculateLocalScore(_stateSets[nodeId], _local_scores[nodeId], _stateSets[index1], _stateSets[index2], _local_scores[index1], _local_scores[index2]);
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
        
        const std::vector<std::string> FlexibleParsimony::getNames()const{
            return _sites.getSequencesNames();
        }
    }
}

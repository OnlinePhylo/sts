#ifndef __hmsbeagle__
#define __hmsbeagle__

#include <string>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include "libhmsbeagle/beagle.h"
#include "phylofunc.hh"
#include <memory>
#include <vector>
#include <stack>
#include <unordered_map>

std::vector< double > get_partials(const std::string& sequence);

class OnlineCalculator
{
public:
    OnlineCalculator() : stateCount(4), instance(-1), next_id(0) {};
    ~OnlineCalculator() {
        beagleFinalizeInstance(instance);
    };

    double calculate_ll( std::shared_ptr< phylo_node > node, std::vector<bool>& visited );
    void initialize(const std::vector<std::string>& seqs);
    int get_id();
    void free_id(int id);
private:
    BeagleInstanceDetails instDetails;
    int stateCount;
    int nPatterns;
    int nPartBuffs;
    int instance;

    int next_id;
    std::stack<int> free_ids;
    std::unordered_map< int, double > id_ll; // caches the root ll at each ID

    int init_beagle();
    void set_eigen_and_rates_and_weights( int instance );
    void grow();
};

int OnlineCalculator::get_id()
{
    if(free_ids.size()>0) {
        int id = free_ids.top();
        free_ids.pop();
        return id;
    }
    if(next_id == nPartBuffs) {
        // oh no! we ran out of buffer slots! need to reallocate
        grow();
    }
    return next_id++;
}

void OnlineCalculator::free_id(int id)
{
    free_ids.push(id);
    id_ll.erase(id);
}

void OnlineCalculator::grow()
{
    nPartBuffs = nPartBuffs * 1.61803398875; // Growing with the golden ratio.
    int new_instance = init_beagle();
    std::cerr << "Growing to " << nPartBuffs << std::endl;

    // Copy all partial probability data into the new instance.
    // XXX: there does not seem to be any way to copy scale factors
    // or to set partials with a particular weight so we need to force a full re-peel.
    double temp[nPatterns * stateCount];
    for(int i=0; i<next_id; i++){
        beagleGetPartials(instance, i, 0, temp);
        beagleSetPartials(new_instance, i, temp);
    }

    // Clear the likelihood cache.
    id_ll.clear();

    set_eigen_and_rates_and_weights( new_instance );    

    // free the old instance
    beagleFinalizeInstance(instance);
    instance = new_instance;
}

///Initialize an instance of the BEAGLE library.

///  \param seqs The sequences.
void OnlineCalculator::initialize(const std::vector<std::string>& seqs)
{

    // Create an instance of the BEAGLE library.
    nPatterns = seqs[1].size();
    nPartBuffs = seqs.size() * 100;
    instance = init_beagle();

    // Add any new sequences to the BEAGLE instance.
    for(int i=0; i<seqs.size(); i++) {
        std::vector< double > seq_partials = get_partials(seqs[i]);
        beagleSetPartials(instance, i, seq_partials.data());
    }
    next_id = seqs.size();

    set_eigen_and_rates_and_weights( instance );
}

int OnlineCalculator::init_beagle()
{
    int new_instance = beagleCreateInstance(
                   0,           /**< Number of tip data elements (input) */
                   nPartBuffs,  /**< Number of partials buffers to create (input) */
                   0,           /**< Number of compact state representation buffers to create (input) */
                   stateCount,  /**< Number of states in the continuous-time Markov chain (input) */
                   nPatterns,   /**< Number of site patterns to be handled by the instance (input) */
                   nPartBuffs,  /**< Number of rate matrix eigen-decomposition buffers to allocate (input) */
                   nPartBuffs,  /**< Number of rate matrix buffers (input) */
                   1,           /**< Number of rate categories (input) */
                   nPartBuffs,  /**< Number of scaling buffers */
                   NULL,        /**< List of potential resource on which this instance is allowed (input, NULL implies no restriction */
                   0,           /**< Length of resourceList list (input) */
                   BEAGLE_FLAG_VECTOR_SSE | BEAGLE_FLAG_PRECISION_DOUBLE | BEAGLE_FLAG_SCALING_AUTO, /**< Bit-flags indicating preferred implementation charactertistics, see BeagleFlags (input) */
                   0,           /**< Bit-flags indicating required implementation characteristics, see BeagleFlags (input) */
                   &instDetails);
    if (new_instance < 0) {
        fprintf(stderr, "Failed to obtain BEAGLE instance\n\n");
        exit(1);
    }
    return new_instance;
}

void OnlineCalculator::set_eigen_and_rates_and_weights( int inst )
{
    // create base frequency array
    double freqs[16] = { 0.25, 0.25, 0.25, 0.25,
                         0.25, 0.25, 0.25, 0.25,
                         0.25, 0.25, 0.25, 0.25,
                         0.25, 0.25, 0.25, 0.25
                       };

    // an eigen decomposition for the JC69 model
    double evec[4 * 4] = {
        1.0,  2.0,  0.0,  0.5,
        1.0,  -2.0,  0.5,  0.0,
        1.0,  2.0, 0.0,  -0.5,
        1.0,  -2.0,  -0.5,  0.0
    };

    double ivec[4 * 4] = {
        0.25,  0.25,  0.25,  0.25,
        0.125,  -0.125,  0.125,  -0.125,
        0.0,  1.0,  0.0,  -1.0,
        1.0,  0.0,  -1.0,  0.0
    };

    double eval[4] = { 0.0, -1.3333333333333333, -1.3333333333333333, -1.3333333333333333 };

    // set the Eigen decomposition
    beagleSetEigenDecomposition(inst, 0, evec, ivec, eval);

    beagleSetStateFrequencies(inst, 0, freqs);
    double rates[1] = { 1 };
    beagleSetCategoryRates(inst, &rates[0]);
    double weights[1] = { 1 };
    beagleSetCategoryWeights(inst, 0, weights);

    double patternWeights[nPatterns];
    for (int i = 0; i < nPatterns; i++) {
        patternWeights[i] = 1.0;
    }
    beagleSetPatternWeights(inst, patternWeights);
}

double OnlineCalculator::calculate_ll( std::shared_ptr< phylo_node > node, std::vector<bool>& visited )
{
    // Prepare the visited vector.
    if(visited.size() < nPartBuffs) {
        visited.resize(nPartBuffs);
    }

    // Accumulate `ops`, a vector of operations, via a depth first search.
    // When likelihoods are cached then operations will only be added for likelihoods that are not cached.
    std::vector< BeagleOperation > ops_tmp, ops;
    std::vector< int > nind; // probability indices
    std::vector< double > lens; // branch lengths
    std::stack< std::shared_ptr< phylo_node > > s;
    s.push(node);
    while(s.size()>0) {
        std::shared_ptr< phylo_node > cur = s.top();
        s.pop();
        if(cur->child1 == NULL){
            visited[cur->id]=true;
            // assert(cur->child2 == NULL);
            continue;
        }
        visited[cur->child1->id]=id_ll.find( cur->child1->id ) != id_ll.end();
        visited[cur->child2->id]=id_ll.find( cur->child2->id ) != id_ll.end();
        s.push(cur->child1);
        s.push(cur->child2);
        if(!visited[cur->id]){
            ops_tmp.push_back( {
                cur->id,           // index of destination, or parent, partials buffer
                BEAGLE_OP_NONE,    // index of scaling buffer to write to (if set to BEAGLE_OP_NONE then calculation of new scalers is disabled)
                BEAGLE_OP_NONE,    // index of scaling buffer to read from (if set to BEAGLE_OP_NONE then use of existing scale factors is disabled)
                cur->child1->id,   // index of first child partials buffer
                cur->child1->id,   // index of transition matrix of first partials child buffer
                cur->child2->id,   // index of second child partials buffer
                cur->child2->id}); // index of transition matrix of second partials child buffer
            nind.push_back(cur->child1->id);
            nind.push_back(cur->child2->id);
            lens.push_back(cur->dist1);
            lens.push_back(cur->dist2);
        }
        visited[cur->id]=true;
    }

    // if we have a cached root LL for this node just return that instead of recalculating
    if( id_ll.find(node->id) != id_ll.end() ){
        return id_ll[ node->id ];
    }

    if(ops_tmp.size() > 0){

        // Need to reverse the order to make post-order.
        ops.insert(ops.begin(), ops_tmp.rbegin(),ops_tmp.rend());

        // Tell BEAGLE to populate the transition matrices for the above edge lengths.
        beagleUpdateTransitionMatrices(instance,     // instance
                                       0,             // eigenIndex
                                       nind.data(),   // probabilityIndices
                                       NULL,          // firstDerivativeIndices
                                       NULL,          // secondDerivativeIndices
                                       lens.data(),   // edgeLengths
                                       nind.size());  // count

        // Create a list of partial likelihood update operations.
        // The order is [dest, destScaling, sourceScaling, source1, matrix1, source2, matrix2].
        // TODO: make this peel only where the new node was added.
        BeagleOperation operations[ops.size()];
        int scaleIndices[ops.size()];
        for(int i=0; i<ops.size(); i++) {
            scaleIndices[i]=ops[i].destinationPartials;
            operations[i] = ops[i];
        }

        // Update the partials.
        beagleUpdatePartials(instance, operations, ops.size(), node->id); // cumulative scaling index
        beagleAccumulateScaleFactors(instance, scaleIndices, ops.size(), BEAGLE_OP_NONE);

    }

    double logL = 0.0;
    int returnCode = 0;

    // Calculate the site likelihoods at the root node.
    int rootIndices[ 1 ] = { node->id };
    int categoryWeightsIndices[ 1 ] = { 0 };
    int stateFrequencyIndices[ 1 ] = { 0 };
    int cumulativeScalingIndices[ 1 ] = { BEAGLE_OP_NONE };
    returnCode = beagleCalculateRootLogLikelihoods(instance, // instance
                     (const int *)rootIndices,               // bufferIndices
                     (const int *)categoryWeightsIndices,    // weights
                     (const int *)stateFrequencyIndices,     // stateFrequencies
                     cumulativeScalingIndices,               // cumulative scaling index
                     1,                                      // count
                     &logL);                                 // outLogLikelihoods

    id_ll[ node->id ] = logL; // stash the log likelihood for later use
    return logL;
}


void add_char( int* table, const int* vals, int len ){
    for(int i=0; i<len; i++)
        table[i] = vals[i];
}

std::vector<double> get_partials(const std::string& sequence)
{
    int n = sequence.size();
    std::vector< double > partials(n * 4);
    const char A = 1 << 0;
    const char C = 1 << 1;
    const char G = 1 << 2;
    const char T = 1 << 3;
    char dna_table[256];
    memset( dna_table, A|C|G|T, 256 );
    dna_table['A']=A;
    dna_table['C']=C;
    dna_table['G']=G;
    dna_table['T']=T;
    dna_table['M']=A|C;
    dna_table['R']=A|G;
    dna_table['W']=A|T;
    dna_table['S']=C|G;
    dna_table['Y']=C|T;
    dna_table['K']=G|T;
    dna_table['V']=A|C|G;
    dna_table['H']=A|C|T;
    dna_table['D']=A|G|T;
    dna_table['B']=C|G|T;
    dna_table['N']=A|C|G|T;
    dna_table['U']=T;
    int k = 0;
    for (int i = 0; i < n; i++) {
        char c = dna_table[ sequence[i] ];
        for(int j=0; j<4; j++){
            partials[k++] = (double)(c & 0x1);
            c >>= 1;
        }
    }
    return partials;
}



#endif //  __hmsbeagle__


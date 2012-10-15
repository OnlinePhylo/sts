#include <string>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <memory>
#include <vector>
#include <stack>
#include <unordered_map>
#include <assert.h>

#include <Bpp/Phyl/Model/AbstractSubstitutionModel.h>
#include <Bpp/Phyl/Model/HKY85.h>
#include <Bpp/Seq/Alphabet/DNA.h>

#include "bpp_shim.hh"
#include "hmsbeagle.hh"

/// Get the ID of an avaiable partial buffer.
/// Allocate more if needed.
///  \return The id.

int OnlineCalculator::get_id()
{
    if(free_ids.size() > 0) {
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

/// Free a partial buffer.
///  \param id The id of the buffer to free.
void OnlineCalculator::free_id(int id)
{
    free_ids.push(id);
    id_ll.erase(id);
}


/// Grow the number of buffers for the BEAGLE instance in OnlineCalculator.
void OnlineCalculator::grow()
{
    nPartBuffs = nPartBuffs * 1.61803398875; // Grow by the golden ratio, to 11 significant figures.
    int new_instance = create_beagle_instance();
    std::cerr << "Growing to " << nPartBuffs << std::endl;

    // Copy all partial probability data into the new instance.
    // XXX: there does not seem to be any way to copy scale factors
    // or to set partials with a particular weight so we need to force a full re-peel.
    double temp[nPatterns * stateCount];
    for(int i = 0; i < next_id; i++) {
        beagleGetPartials(instance, i, 0, temp);
        beagleSetPartials(new_instance, i, temp);
    }

    // Clear the likelihood cache.
    id_ll.clear();

    set_eigen_and_rates_and_weights(new_instance);

    // Free the old instance.
    beagleFinalizeInstance(instance);
    instance = new_instance;
}

///Initialize an instance of the BEAGLE library with partials coming from sequences.
///  \param seqs A string vector of sequences.
void OnlineCalculator::initialize(const std::vector<std::string>& seqs)
{
    // Create an instance of the BEAGLE library.
    assert(!seqs.empty());
    nPatterns = seqs[0].size();
    nPartBuffs = seqs.size() * 100; // AD: wouldn't it make sense to make this an argument?
    instance = create_beagle_instance();

    // Add the sequences to the BEAGLE instance.
    for(int i = 0; i < seqs.size(); i++) {
        std::vector< double > seq_partials = get_partials(seqs[i]);
        beagleSetPartials(instance, i, seq_partials.data());
    }
    next_id = seqs.size();

    // XXX Hard coding a model for the time being. Note that grow still uses old call.
    bpp::DNA alphabet;
    bpp::HKY85 model(&alphabet);
    set_eigen_and_rates_and_weights(instance, model);
}

/// Create an instance of the BEAGLE library.
/// \return An integer identifier for the instance.
int OnlineCalculator::create_beagle_instance()
{
    int new_instance =
        beagleCreateInstance(
            0,           // Number of tip data elements (input)
            nPartBuffs,  // Number of partials buffers to create (input)
            0,           // Number of compact state representation buffers to create (input)
            stateCount,  // Number of states in the continuous-time Markov chain (input)
            nPatterns,   // Number of site patterns to be handled by the instance (input)
            nPartBuffs,  // Number of rate matrix eigen-decomposition buffers to allocate (input)
            nPartBuffs,  // Number of rate matrix buffers (input)
            1,           // Number of rate categories (input)
            nPartBuffs,  // Number of scaling buffers
            NULL,        // List of potential resource on which this instance is allowed (input, NULL implies no restriction
            0,           // Length of resourceList list (input)
            BEAGLE_FLAG_VECTOR_SSE | BEAGLE_FLAG_PRECISION_DOUBLE | BEAGLE_FLAG_SCALING_AUTO,
            // Bit-flags indicating preferred implementation charactertistics, see BeagleFlags (input)
            0,           // Bit-flags indicating required implementation characteristics, see BeagleFlags (input)
            &instDetails);
    if(new_instance < 0) {
        fprintf(stderr, "Fatal: failed to obtain BEAGLE instance.\n\n");
        exit(1);
    }
    return new_instance;
}

void OnlineCalculator::set_eigen_and_rates_and_weights(int inst)
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
    for(int i = 0; i < nPatterns; i++) {
        patternWeights[i] = 1.0;
    }
    beagleSetPatternWeights(inst, patternWeights);
}

void
OnlineCalculator::set_eigen_and_rates_and_weights(int inst, const bpp::AbstractReversibleSubstitutionModel& model)
{
    double freqs[4], evec[16], ivec[16], eval[4]; // XXX Fix for transition to AA
    blit_matrix_to_array(ivec, model.getRowLeftEigenVectors()); // inverse eigenvectors
    blit_matrix_to_array(evec, model.getColumnRightEigenVectors()); // eigenvectors
    blit_vector_to_array(eval, model.getEigenValues());
    beagleSetEigenDecomposition(inst, 0, evec, ivec, eval);

    blit_vector_to_array(freqs, model.getFrequencies());
    beagleSetStateFrequencies(inst, 0, freqs);

    double rate = model.getRate(), weight = 1;
    beagleSetCategoryRates(inst, &rate);
    beagleSetCategoryWeights(inst, 0, &weight);
}

/// Calculate the log likelihood
/// \param node The root std::shared_ptr<phylo_node> at which to start computation.
/// \param visited A std::vector<bool>& with enough entries to store the visited status of all daughter nodes.
/// \return the log likelihood.
double OnlineCalculator::calculate_ll(std::shared_ptr< phylo_node > node, std::vector<bool>& visited)
{
    // Resize if visited vector is not big enough.
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
    // Recursively traverse the tree, accumulating operations.
    while(s.size() > 0) {
        std::shared_ptr< phylo_node > cur = s.top(); // AD: seems like we could use a weak_ptr here, right? Wouldn't that be imperceptibly faster?
        s.pop();
        if(cur->is_leaf()) {
            // We are at a leaf.
            assert(cur->child2->node == NULL);
            visited[cur->id] = true;
            continue;
        }
        // Mark a child as visited if we have already calculated its log likelihood.
        visited[cur->child1->node->id] = id_ll.find(cur->child1->node->id) != id_ll.end();
        visited[cur->child2->node->id] = id_ll.find(cur->child2->node->id) != id_ll.end();
        // AD: do we not assume that children of visited nodes have themselves been visited? Seems like we could avoid
        // these pushes if so.
        s.push(cur->child1->node);
        s.push(cur->child2->node);
        if(!visited[cur->id]) {
            ops_tmp.push_back(BeagleOperation({
                cur->id,           // index of destination, or parent, partials buffer
                BEAGLE_OP_NONE,    // index of scaling buffer to write to (if set to BEAGLE_OP_NONE then calculation of new scalers is disabled)
                BEAGLE_OP_NONE,    // index of scaling buffer to read from (if set to BEAGLE_OP_NONE then use of existing scale factors is disabled)
                cur->child1->node->id,   // index of first child partials buffer
                cur->child1->node->id,   // index of transition matrix of first partials child buffer
                cur->child2->node->id,   // index of second child partials buffer
                cur->child2->node->id    // index of transition matrix of second partials child buffer
            }));
            nind.push_back(cur->child1->node->id);
            nind.push_back(cur->child2->node->id);
            lens.push_back(cur->child1->length);
            lens.push_back(cur->child2->length);
        }
        visited[cur->id] = true;
    }

    // If we have a cached root LL for this node just return that instead of recalculating.
    if(id_ll.find(node->id) != id_ll.end()) {
        return id_ll[ node->id ];
    }

    if(ops_tmp.size() > 0) { // If we actually need to do some operations.

        // Reverse the order of operations to make them post-order.
        ops.insert(ops.begin(), ops_tmp.rbegin(), ops_tmp.rend());

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
        for(int i = 0; i < ops.size(); i++) {
            scaleIndices[i] = ops[i].destinationPartials;
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
                 &logL);                                 // OUT: log likelihood

    id_ll[ node->id ] = logL; // Record the log likelihood for later use.
    return logL;
}

// Turn a sequence into partials.
std::vector<double> get_partials(const std::string& sequence)
{
    int n = sequence.size();
    std::vector< double > partials(n * 4);
    const char A = 1 << 0;
    const char C = 1 << 1;
    const char G = 1 << 2;
    const char T = 1 << 3;
    char dna_table[256];
    memset(dna_table, A | C | G | T, 256);
    dna_table['A'] = A;
    dna_table['C'] = C;
    dna_table['G'] = G;
    dna_table['T'] = T;
    dna_table['M'] = A | C;
    dna_table['R'] = A | G;
    dna_table['W'] = A | T;
    dna_table['S'] = C | G;
    dna_table['Y'] = C | T;
    dna_table['K'] = G | T;
    dna_table['V'] = A | C | G;
    dna_table['H'] = A | C | T;
    dna_table['D'] = A | G | T;
    dna_table['B'] = C | G | T;
    dna_table['N'] = A | C | G | T;
    dna_table['U'] = T;
    int k = 0;
    for(int i = 0; i < n; i++) {
        char c = dna_table[ toupper(sequence[i]) ];
        for(int j = 0; j < 4; j++) {
            partials[k++] = (double)(c & 0x1);
            c >>= 1;
        }
    }
    return partials;
}


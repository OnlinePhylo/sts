#include <string>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <memory>
#include <vector>
#include <stack>
#include <unordered_map>
#include <assert.h>
#include <iostream>

#include <Bpp/Phyl/Model/JCnuc.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Seq/Alphabet/DNA.h>

#include "bpp_shim.hh"
#include "hmsbeagle.hh"

namespace sts
{

/// Get the ID of an avaiable partial buffer.
/// Allocate more if needed.
///  \return The id.

int online_calculator::get_id()
{
    if(free_ids.size() > 0) {
        int id = free_ids.top();
        free_ids.pop();
        return id;
    }
    if(next_id == num_buffers) {
        // oh no! we ran out of buffer slots! need to reallocate
        grow();
    }
    return next_id++;
}

/// Free a partial buffer.
///  \param id The id of the buffer to free.
void online_calculator::free_id(int id)
{
    free_ids.push(id);
    map_id_ll.erase(id);
}


/// Grow the number of buffers for the BEAGLE instance in online_calculator.
void online_calculator::grow()
{
    num_buffers = num_buffers * 1.61803398875; // Grow by the golden ratio, to 11 significant figures.
    int new_instance = create_beagle_instance();

    // Copy all partial probability data into the new instance.
    // There does not seem to be any way to copy scale factors
    // or to set partials with a particular weight index so we need to force a full re-peel.
    double temp[sites->getNumberOfSites() * sites->getAlphabet()->getSize()];
    for(int i = 0; i < next_id; i++) {
        beagleGetPartials(instance, i, 0, temp);
        beagleSetPartials(new_instance, i, temp);
    }

    // Clear the likelihood cache.
    map_id_ll.clear();
    set_eigen_and_rates_and_weights(new_instance);

    // Free the old instance.
    beagleFinalizeInstance(instance);
    instance = new_instance;
}

///Initialize an instance of the BEAGLE library with partials coming from sequences.
///  \param sites The sites
///  \param model Substitution model
void online_calculator::initialize(std::shared_ptr<bpp::SiteContainer> sites, std::shared_ptr<bpp::SubstitutionModel> model)
{
    this->sites = sites;
    this->model = model;

    const int n_seqs = sites->getNumberOfSequences();

    // Create an instance of the BEAGLE library.
    assert(sites->getNumberOfSites());
    assert(n_seqs);

    num_buffers = n_seqs * 100; // AD: wouldn't it make sense to make this an argument?

    instance = create_beagle_instance();

    // Add the sequences to the BEAGLE instance.
    for(int i = 0; i < n_seqs; i++) {
        std::vector<double> seq_partials = get_partials(sites->getSequence(i), *model, sites->getAlphabet());
        beagleSetPartials(instance, i, seq_partials.data());
    }
    next_id = n_seqs;

    set_eigen_and_rates_and_weights(instance);

    initialized = true;
}


/// Create an instance of the BEAGLE library.
/// \return An integer identifier for the instance.
int online_calculator::create_beagle_instance()
{
    int new_instance =
        beagleCreateInstance(
            0,           // Number of tip data elements (input)
            num_buffers,  // Number of partials buffers to create (input)
            0,           // Number of compact state representation buffers to create (input)
            sites->getAlphabet()->getSize(),  // Number of states in the continuous-time Markov chain (input)
            sites->getNumberOfSites(),   // Number of site patterns to be handled by the instance (input)
            num_buffers,  // Number of rate matrix eigen-decomposition buffers to allocate (input)
            num_buffers,  // Number of rate matrix buffers (input)
            1,           // Number of rate categories (input)
            num_buffers,  // Number of scaling buffers
            NULL,        // List of potential resource on which this instance is allowed (input, NULL implies no restriction
            0,           // Length of resourceList list (input)
            BEAGLE_FLAG_VECTOR_SSE | BEAGLE_FLAG_PRECISION_DOUBLE | BEAGLE_FLAG_SCALING_AUTO,
            // Bit-flags indicating preferred implementation charactertistics, see BeagleFlags (input)
            0,           // Bit-flags indicating required implementation characteristics, see BeagleFlags (input)
            &instance_details);
    if(new_instance < 0) {
        fprintf(stderr, "Fatal: failed to obtain BEAGLE instance.\n\n");
        exit(1);
    }
    return new_instance;
}

/// Set eigen rates and weights
///   \param inst Beagle instance
void online_calculator::set_eigen_and_rates_and_weights(int inst)
{
    assert(model);
    int n_states = model->getAlphabet()->getSize();
    int n_patterns = sites->getNumberOfSites();
    double evec[n_states * n_states], ivec[n_states * n_states], eval[n_states];

    blit_matrix_to_array(ivec, model->getRowLeftEigenVectors()); // inverse eigenvectors
    blit_matrix_to_array(evec, model->getColumnRightEigenVectors()); // eigenvectors
    blit_vector_to_array(eval, model->getEigenValues());
    beagleSetEigenDecomposition(inst, 0, evec, ivec, eval);

    beagleSetStateFrequencies(inst, 0, model->getFrequencies().data());

    double rate = model->getRate(), weight = 1;
    beagleSetCategoryRates(inst, &rate);
    beagleSetCategoryWeights(inst, 0, &weight);

    double patternWeights[n_patterns];
    for(int i = 0; i < n_patterns; i++) {
        patternWeights[i] = 1.0;
    }
    beagleSetPatternWeights(inst, patternWeights);
}

/// Calculate the log likelihood
/// \param node The root std::shared_ptr<phylo_node> at which to start computation.
/// \param visited A std::vector<bool>& with enough entries to store the visited status of all daughter nodes.
/// \return the log likelihood.
double online_calculator::calculate_ll(std::shared_ptr< phylo_node > node, std::vector<bool>& visited)
{
    // Resize if visited vector is not big enough.
    if(visited.size() < num_buffers) {
        visited.resize(num_buffers);
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
            visited[cur->id] = true;
            continue;
        }
        // Mark a child as visited if we have already calculated its log likelihood.
        visited[cur->child1->node->id] = map_id_ll.find(cur->child1->node->id) != map_id_ll.end();
        visited[cur->child2->node->id] = map_id_ll.find(cur->child2->node->id) != map_id_ll.end();
        // AD: do we not assume that children of visited nodes have themselves been visited? Seems like we could avoid
        // these pushes if so.
        // Reply: Keeping these in supports full peeling when needed, e.g. after reallocating the beagle instance because
        // we needed to grow.
        s.push(cur->child1->node);
        s.push(cur->child2->node);
        if(!visited[cur->id]) {
            ops_tmp.push_back(BeagleOperation( {
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
    if(map_id_ll.find(node->id) != map_id_ll.end()) {
        return map_id_ll[ node->id ];
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

    map_id_ll[ node->id ] = logL; // Record the log likelihood for later use.
    return logL;
}

} // namespace sts

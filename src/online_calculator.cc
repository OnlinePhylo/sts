/// \file online_calculator.cc
/// \author metatangle, inc.
/// \brief Interface between STS and beagle-lib

#include "online_calculator.h"

#include <cassert>
#include <iostream>

#include "bpp_shim.h"
#include "node.h"
#include "node_ptr.h"
#include "edge.h"

namespace sts
{
namespace likelihood
{

std::string beagle_errstring(const int beagle_error_code)
{
    switch(beagle_error_code) {
        case BEAGLE_SUCCESS:                      return "BEAGLE_SUCCESS";
        case BEAGLE_ERROR_GENERAL:                return "BEAGLE_ERROR_GENERAL";
        case BEAGLE_ERROR_OUT_OF_MEMORY:          return "BEAGLE_ERROR_OUT_OF_MEMORY";
        case BEAGLE_ERROR_UNIDENTIFIED_EXCEPTION: return "BEAGLE_ERROR_UNIDENTIFIED_EXCEPTION";
        case BEAGLE_ERROR_UNINITIALIZED_INSTANCE: return "BEAGLE_ERROR_UNINITIALIZED_INSTANCE";
        case BEAGLE_ERROR_OUT_OF_RANGE:           return "BEAGLE_ERROR_OUT_OF_RANGE";
        case BEAGLE_ERROR_NO_RESOURCE:            return "BEAGLE_ERROR_NO_RESOURCE";
        case BEAGLE_ERROR_NO_IMPLEMENTATION:      return "BEAGLE_ERROR_NO_IMPLEMENTATION";
        case BEAGLE_ERROR_FLOATING_POINT:         return "BEAGLE_ERROR_FLOATING_POINT";
        default: return "Unknown Beagle error: " + std::to_string(beagle_error_code);
    }
}

Online_calculator::~Online_calculator()
{
    if(instance >= 0)
        beagleFinalizeInstance(instance);
}

/// Get the ID of an available partial buffer.
/// Allocate more if needed.
///  \return The id.

int Online_calculator::get_id()
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
void Online_calculator::free_id(int id)
{
    free_ids.push(id);
}


/// Grow the number of buffers for the BEAGLE instance in Online_calculator.
void Online_calculator::grow()
{
    num_buffers = num_buffers * 1.61803398875; // Grow by the golden ratio, to 11 significant figures.
    int new_instance = create_beagle_instance();

    // Copy all partial probability data into the new instance.
    // There does not seem to be any way to copy scale factors
    // or to set partials with a particular weight index so we need to force a full re-peel.
    const std::unique_ptr<double[]> temp(new double[sites->getNumberOfSites() * sites->getAlphabet()->getSize()]);
    for(int i = 0; i < next_id; i++) {
        beagleGetPartials(instance, i, 0, temp.get());
        beagleSetPartials(new_instance, i, temp.get());
    }

    // Clear the likelihood cache.
    node_ll_map.clear();
    set_eigen_and_rates_and_weights(new_instance);

    // Free the old instance.
    beagleFinalizeInstance(instance);
    instance = new_instance;
}

void Online_calculator::register_node(sts::particle::Node_ptr n)
{
    assert(node_buffer_map.count(n.get()) == 0);
    node_buffer_map[n.get()] = get_id();
}

/// Register a leaf node with the calculator

/// \param n Node
/// \param taxon The taxon name for a node. Must match one of the input taxa passed to
///              sts::likelihood::Online_calculator::initialize.
void Online_calculator::register_leaf(sts::particle::Node_ptr n, const std::string taxon)
{
    assert(node_buffer_map.count(n.get()) == 0);
    assert(taxon_buffer_map.count(taxon) == 1);
    assert(n->is_leaf());
    node_buffer_map[n.get()] = taxon_buffer_map[taxon];
}

int Online_calculator::get_buffer(sts::particle::Node_ptr n)
{
    if(node_buffer_map.count(n.get()) == 0) register_node(n);
    return node_buffer_map[n.get()];
}

void Online_calculator::unregister_node(const sts::particle::Node* n)
{
    if(node_buffer_map.count(n) == 0) return;
    free_id(node_buffer_map[n]);
    node_buffer_map.erase(n);
    node_ll_map.erase(n);
}

///Initialize an instance of the BEAGLE library with partials coming from sequences.
///  \param sites The sites
///  \param model Substitution model
void Online_calculator::initialize(std::shared_ptr<bpp::SiteContainer> sites, std::shared_ptr<bpp::SubstitutionModel> model)
{
    this->sites = sites;
    this->model = model;
    this->weights = std::vector<double>(sites->getNumberOfSites(), 1.0);

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
        taxon_buffer_map[sites->getSequencesNames()[i]] = i;
    }
    next_id = n_seqs;

    set_eigen_and_rates_and_weights(instance);
}


/// Create an instance of the BEAGLE library.
/// \return An integer identifier for the instance.
int Online_calculator::create_beagle_instance()
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
            num_buffers + 1, // Number of scaling buffers
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
void Online_calculator::set_eigen_and_rates_and_weights(int inst)
{
    assert(model);
    int n_states = model->getAlphabet()->getSize();
    const std::unique_ptr<double[]>
    evec(new double[n_states * n_states]),
         ivec(new double[n_states * n_states]),
         eval(new double[n_states]);

    blit_matrix_to_array(ivec.get(), model->getRowLeftEigenVectors()); // inverse eigenvectors
    blit_matrix_to_array(evec.get(), model->getColumnRightEigenVectors()); // eigenvectors
    blit_vector_to_array(eval.get(), model->getEigenValues());
    beagleSetEigenDecomposition(inst, 0, evec.get(), ivec.get(), eval.get());

    beagleSetStateFrequencies(inst, 0, model->getFrequencies().data());

    double rate = model->getRate(), weight = 1;
    beagleSetCategoryRates(inst, &rate);
    beagleSetCategoryWeights(inst, 0, &weight);

    beagleSetPatternWeights(inst, weights.data());
}

/// Set the weights
/// \param weights A vector with length equal to the number of sites
void Online_calculator::set_weights(std::vector<double> weights)
{
    this->weights = weights;
    assert(model);
    assert(weights.size() == sites->getNumberOfSites());

    beagleSetPatternWeights(instance, weights.data());
}

/// Calculate the log likelihood
/// \param node The root std::shared_ptr<sts::particle::Node> at which to start computation.
/// \param visited A std::vector<bool>& with enough entries to store the visited status of all daughter nodes.
/// \return the log likelihood.
double Online_calculator::calculate_ll(sts::particle::Node_ptr node, std::unordered_set<sts::particle::Node_ptr>& visited)
{
    // Accumulate `ops`, a vector of operations, via a depth first search.
    // When likelihoods are cached then operations will only be added for likelihoods that are not cached.
    std::vector<BeagleOperation> ops_tmp, ops;
    std::vector<int> nind; // probability indices
    std::vector<double> lens; // branch lengths
    std::stack<sts::particle::Node_ptr> s;
    s.push(node);
    // Recursively traverse the tree, accumulating operations.
    while(s.size() > 0) {
        sts::particle::Node_ptr cur = s.top();
        s.pop();
        if(cur->is_leaf()) {
            // We are at a leaf.
            visited.insert(cur);
            continue;
        }
        // Mark a child as visited if we have already calculated its log likelihood.
        if(node_ll_map.count(cur->child1->node.get()) != 0) visited.insert(cur->child1->node);
        if(node_ll_map.count(cur->child2->node.get()) != 0) visited.insert(cur->child2->node);
        // AD: do we not assume that children of visited nodes have themselves been visited? Seems like we could avoid
        // these pushes if so.
        // Reply: Keeping these in supports full peeling when needed, e.g. after reallocating the beagle instance because
        // we needed to grow.
        s.push(cur->child1->node);
        s.push(cur->child2->node);
        if(visited.count(cur) == 0) {
            ops_tmp.push_back(BeagleOperation( {
                get_buffer(cur),           // index of destination, or parent, partials buffer
                BEAGLE_OP_NONE,    // index of scaling buffer to write to (if set to BEAGLE_OP_NONE then calculation of new scalers is disabled)
                BEAGLE_OP_NONE,    // index of scaling buffer to read from (if set to BEAGLE_OP_NONE then use of existing scale factors is disabled)
                get_buffer(cur->child1->node),   // index of first child partials buffer
                get_buffer(cur->child1->node),   // index of transition matrix of first partials child buffer
                get_buffer(cur->child2->node),   // index of second child partials buffer
                get_buffer(cur->child2->node)    // index of transition matrix of second partials child buffer
            }));
            nind.push_back(get_buffer(cur->child1->node));
            nind.push_back(get_buffer(cur->child2->node));
            lens.push_back(cur->child1->length);
            lens.push_back(cur->child2->length);
        }
        visited.insert(cur);
    }

    // If we have a cached root LL for this node just return that instead of recalculating.
    if(!verify_cached_ll && node_ll_map.count(node.get()) != 0) {
        return node_ll_map[node.get()];
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
        const std::unique_ptr<BeagleOperation[]> operations(new BeagleOperation[ops.size()]);
        const std::unique_ptr<int[]> scaleIndices(new int[ops.size()]);
        for(size_t i = 0; i < ops.size(); i++) {
            scaleIndices[i] = ops[i].destinationPartials;
            operations[i] = ops[i];
        }

        // Update the partials.
        beagleUpdatePartials(instance, operations.get(), ops.size(), get_buffer(node)); // cumulative scaling index
        beagleAccumulateScaleFactors(instance, scaleIndices.get(), ops.size(), num_buffers);

    }

    double logL = 0.0;
    int returnCode = 0;

    // Calculate the site likelihoods at the root node.
    int rootIndices[ 1 ] = { get_buffer(node) };
    int categoryWeightsIndices[ 1 ] = { 0 };
    int stateFrequencyIndices[ 1 ] = { 0 };
    int cumulativeScalingIndices[ 1 ] = { num_buffers };
    returnCode = beagleCalculateRootLogLikelihoods(instance, // instance
                 (const int *)rootIndices,               // bufferIndices
                 (const int *)categoryWeightsIndices,    // weights
                 (const int *)stateFrequencyIndices,     // stateFrequencies
                 cumulativeScalingIndices,               // cumulative scaling index
                 1,                                      // count
                 &logL);                                 // OUT: log likelihood

    if(returnCode)
        throw std::runtime_error(beagle_errstring(returnCode));
    // Verify LL if requested.
    if(verify_cached_ll && node_ll_map.count(node.get()))
        assert(std::abs(node_ll_map[node.get()] - logL) < 1e-5);

    node_ll_map[node.get()] = logL; // Record the log likelihood for later use.
    return logL;
}

/// Invalidate the cached log likelihood for a node.

/// \param n The node to invalidate.
void Online_calculator::invalidate(std::shared_ptr< sts::particle::Node > n)
{
    node_ll_map.erase(n.get());
}

} // namespace likelihood
} // namespace sts

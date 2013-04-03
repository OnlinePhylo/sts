/// \file beagle_tree_likelihood.cc

#include "beagle_tree_likelihood.h"
#include "bpp_shim.h"
#include "likelihood_vector.h"
#include "util.h"
#include "online_util.h"

#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Seq/Container/SiteContainer.h>

#include <cassert>
#include <stdexcept>
#include <stack>
#include <forward_list>

#include <iostream>

#include <unordered_set>


using sts::likelihood::blit_vector_to_array;
using sts::likelihood::blit_matrix_to_array;
using sts::likelihood::get_partials;
using sts::util::beagle_check;

namespace sts { namespace online {

template<typename T> void hash_combine(size_t& seed, const T& v)
{
    std::hash<T> h;
    seed ^= h(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}


/// Hash a node pointer using a combination of the address,
/// sons addresses, and distance to father
size_t hash_node(const bpp::Node* node)
{
    std::hash<const bpp::Node*> h;
    size_t seed = h(node);
    if(node->hasDistanceToFather())
        hash_combine(seed, node->getDistanceToFather());
    for(size_t i = 0; i < node->getNumberOfSons(); i++)
        hash_combine(seed, node->getSon(i));
    return seed;
}

/// Mark nodes and **ascendants** dirty if changed
template<typename N>
std::vector<const N*> postorder_find_changed(const bpp::TreeTemplate<N>& tree, const std::unordered_map<const N*, size_t>& state)
{
    std::unordered_set<const N*> dirty;
    std::vector<const N*> postorder_nodes = postorder(tree.getRootNode());
    for(const N* node : postorder_nodes) {
        auto it = state.find(node);
        bool is_dirty = it == state.end() ||
                        it->second != hash_node(node) ||
                        dirty.find(node) != dirty.end();
        if(is_dirty) {
            dirty.insert(node);
            dirty.insert(node->getFather());
        }
    }
    auto it = std::remove_if(postorder_nodes.begin(), postorder_nodes.end(),
                             [&dirty](const N* node) { return dirty.find(node) == dirty.end(); });
    postorder_nodes.erase(it, postorder_nodes.end());
    return postorder_nodes;
}

/// Mark nodes and **descendants** dirty if changed
template<typename N>
std::vector<const N*> preorder_find_changed(const bpp::TreeTemplate<N>& tree, const std::unordered_map<const N*, size_t>& state)
{
    std::unordered_set<const N*> dirty;
    std::vector<const N*> preorder_nodes = preorder(tree.getRootNode());
    for(const N* node : preorder_nodes) {
        if(node == tree.getRootNode())
            continue;
        auto it = state.find(node);
        bool is_dirty = it == state.end() ||
                        it->second != hash_node(node) ||
                        dirty.find(node) != dirty.end();
        if(is_dirty) {
            dirty.insert(node);
            for(size_t i = 0; i < node->getNumberOfSons(); i++)
                dirty.insert(node->getSon(i));
            // If node is dirty, and father is root, mark sibling
            if(node->getFather() == tree.getRootNode())
                dirty.insert(siblings(node)[0]);
        }
    }
    auto it = std::remove_if(preorder_nodes.begin(), preorder_nodes.end(),
                             [&dirty](const N* node) { return dirty.find(node) == dirty.end(); });
    preorder_nodes.erase(it, preorder_nodes.end());
    return preorder_nodes;
}


BeagleTreeLikelihood::BeagleTreeLikelihood(const bpp::SiteContainer& sites,
                                               const bpp::SubstitutionModel& model,
                                               const bpp::DiscreteDistribution& rate_dist,
                                               const size_t n_scratch_buffers) :
    beagleInstance_(-1),
    nSites_(sites.getNumberOfSites()),
    nStates_(model.getNumberOfStates()),
    nRates_(rate_dist.getNumberOfCategories()),
    nSeqs_(sites.getNumberOfSequences()),
    // Allocate two buffers for each node in the tree (to store distal, proximal vectors)
    // plus `scratch_buffer_count` BONUS buffers
    nBuffers_(4 * nSeqs_ - 2 + n_scratch_buffers),
    nScratchBuffers_(n_scratch_buffers),
    rateDist(&rate_dist),
    model(&model),
    tree(nullptr)
{
    assert(nRates_ >= 1);

    leafBuffer.reserve(nSeqs_);

    beagleInstance_ = beagleCreateInstance(
            0,              // Number of tip data elements (input)
            nBuffers_,      // Number of partials buffers to create (input)
            0,              // Number of compact state representation buffers to create (input)
            nStates_,       // Number of states in the continuous-time Markov chain (input)
            nSites_,        // Number of site patterns to be handled by the instance (input)
            1,              // Number of rate matrix eigen-decomposition buffers to allocate (input)
            nBuffers_,      // Number of rate matrix buffers (input)
            nRates_,        // Number of rate categories (input)
            nBuffers_ + 2,  // Number of scaling buffers - 1 extra buffer for prox, distal
            NULL,           // List of potential resource on which this instance is allowed (input, NULL implies no
                            // restriction
            0,              // Length of resourceList list (input)
            BEAGLE_FLAG_VECTOR_SSE | BEAGLE_FLAG_PRECISION_DOUBLE | BEAGLE_FLAG_SCALING_AUTO, // Bit-flags indicating
                            //preferred implementation charactertistics, see BeagleFlags (input)
            0,              // Bit-flags indicating required implementation characteristics, see BeagleFlags (input)
            &instanceDetails);
    if(beagleInstance_ < 0)
        beagle_check(beagleInstance_);

    // Load tips
    for(size_t i = 0; i < nSeqs_; i++)
        registerLeaf(sites.getSequence(i));

    // Weight all sites equally -
    // for online inference, we don't want to compress sites.
    std::vector<double> pattern_weights(sites.getNumberOfSites(), 1.0);
    beagle_check(beagleSetPatternWeights(beagleInstance_, pattern_weights.data()));
}

BeagleTreeLikelihood::BeagleTreeLikelihood(BeagleTreeLikelihood&& other) :
    beagleInstance_(other.beagleInstance_),
    nSites_(other.nSites_),
    nStates_(other.nStates_),
    nRates_(other.nRates_),
    nSeqs_(other.nSeqs_),
    nBuffers_(other.nBuffers_),
    nScratchBuffers_(other.nScratchBuffers_),
    scratchBuffers_(std::move(other.scratchBuffers_)),
    instanceDetails(std::move(other.instanceDetails)),
    leafBuffer(std::move(other.leafBuffer)),
    rateDist(other.rateDist),
    model(other.model),
    distalNodeBuffer(std::move(other.distalNodeBuffer)),
    proxNodeBuffer(std::move(other.proxNodeBuffer)),
    distalNodeState(std::move(other.distalNodeState)),
    proxNodeState(std::move(other.proxNodeState))
{ }

BeagleTreeLikelihood::~BeagleTreeLikelihood()
{
    if(beagleInstance_ >= 0)
        beagleFinalizeInstance(beagleInstance_);
}

void BeagleTreeLikelihood::initialize(const bpp::SubstitutionModel& model,
                                        const bpp::DiscreteDistribution& rate_dist,
                                        bpp::TreeTemplate<bpp::Node>& tree)
{
    this->rateDist = &rate_dist;
    this->model = &model;
    this->tree = &tree;
    verifyInitialized();

    distalNodeBuffer.clear();
    proxNodeBuffer.clear();
    distalNodeState.clear();
    proxNodeState.clear();

    loadRateDistribution(rate_dist);
    loadSubstitutionModel(model);

    // Slide root to far right side of root branch
    tree.getRootNode()->getSon(0)->setDistanceToFather(tree.getRootNode()->getSon(0)->getDistanceToFather() +
                                                       tree.getRootNode()->getSon(1)->getDistanceToFather());
    tree.getRootNode()->getSon(1)->setDistanceToFather(0.0);

    // Fill buffer maps
    size_t buffer = nSeqs_;
    const std::vector<bpp::Node*> nodes = tree.getNodes();
    // Distal buffer
    for(const bpp::Node* n : nodes) {
        if(n->isLeaf()) {
            const std::string& name = n->getName();
            assert(leafBuffer.count(name) > 0);
            distalNodeBuffer[n] = leafBuffer.at(name);
        } else {
            assert(buffer < nBuffers_);
            assert(n->getNumberOfSons() == 2);
            assert(distalNodeBuffer.find(n) == distalNodeBuffer.end());
            distalNodeBuffer[n] = buffer;
            buffer++;
        }
    }

    // Proximal buffer
    const bpp::Node* root = tree.getRootNode();
    for(size_t i = 0; i < root->getNumberOfSons(); i++) {
        const bpp::Node* n = root->getSon(i);
        proxNodeBuffer[n] = distalNodeBuffer.at(siblings(n)[0]);
    }
    for(const bpp::Node* n : nodes) {
        if(n->isLeaf() || n == root)
            continue;
        assert(buffer < nBuffers_);
        assert(n->getNumberOfSons() == 2);
        for(size_t i = 0; i < 2; ++i) {
            const bpp::Node* son = n->getSon(i);
            assert(distalNodeBuffer.find(son) != distalNodeBuffer.end());
            assert(proxNodeBuffer.find(son) == proxNodeBuffer.end());
            proxNodeBuffer[son] = buffer;
            buffer++;
        }
    }

    // Store indices for some scratch buffers
    scratchBuffers_.clear();
    for(size_t i = 0; i < nScratchBuffers_; i++) {
        assert(buffer < nBuffers_);
        scratchBuffers_.push_back(buffer++);
    }
}

size_t BeagleTreeLikelihood::registerLeaf(const bpp::Sequence& sequence)
{
    verifyInitialized();
    int buffer = leafBuffer.size();
    if(leafBuffer.count(sequence.getName()) > 0)
        throw std::runtime_error("Duplicate sequence name: " + sequence.getName());

    leafBuffer[sequence.getName()] = buffer;
    const std::vector<double> seq_partials = get_partials(sequence, *model, nRates_);
    assert(seq_partials.size() == sequence.size() * nStates_ * nRates_);
    beagle_check(beagleSetPartials(beagleInstance_, buffer, seq_partials.data()));
    return buffer;
}

void BeagleTreeLikelihood::loadSubstitutionModel(const bpp::SubstitutionModel& model)
{
    verifyInitialized();
    // Eigendecomposition
    const std::unique_ptr<double[]> evec(new double[nStates_ * nStates_]),
                                    ivec(new double[nStates_ * nStates_]),
                                    eval(new double[nStates_]);
    blit_matrix_to_array(ivec.get(), model.getRowLeftEigenVectors());     // inverse eigenvectors
    blit_matrix_to_array(evec.get(), model.getColumnRightEigenVectors()); // eigenvectors
    blit_vector_to_array(eval.get(), model.getEigenValues());
    beagle_check(beagleSetEigenDecomposition(beagleInstance_, 0, evec.get(), ivec.get(), eval.get()));
    // State frequencies
    beagle_check(beagleSetStateFrequencies(beagleInstance_, 0, model.getFrequencies().data()));
}

void BeagleTreeLikelihood::loadRateDistribution(const bpp::DiscreteDistribution& rate_dist)
{
    assert(rate_dist.getNumberOfCategories() == nRates_ &&
           "Unexpected rate category count");
    verifyInitialized();
    const std::vector<double>& categories = rate_dist.getCategories();
    const std::vector<double>& weights = rate_dist.getProbabilities();

    beagle_check(beagleSetCategoryRates(beagleInstance_, categories.data()));
    beagle_check(beagleSetCategoryWeights(beagleInstance_, 0, weights.data()));
}


void BeagleTreeLikelihood::calculateDistalPartials()
{
    verifyInitialized();
    std::vector<const bpp::Node*> postorder_nodes = postorder_find_changed(*tree, distalNodeState);

    // For tracking BEAGLE operations
    std::vector<BeagleOperation> operations;
    std::vector<int> node_indices;      // probability indices
    std::vector<double> branch_lengths;

    // Traverse nodes in postorder, adding BeagleOperations to update each
    for(const bpp::Node* n : postorder_nodes) {
        if(n->isLeaf()) {
            const std::string& name = n->getName();
            assert(leafBuffer.count(name) > 0);
            distalNodeBuffer[n] = leafBuffer.at(name);
        } else {
            assert(n->getNumberOfSons() == 2);
            assert(distalNodeBuffer.find(n) != distalNodeBuffer.end());
            for(size_t i = 0; i < 2; ++i) {
                assert(distalNodeBuffer.count(n->getSon(i)) > 0);
            }
            int buffer = distalNodeBuffer.at(n);
            int child1_buffer = distalNodeBuffer.at(n->getSon(0)),
                child2_buffer = distalNodeBuffer.at(n->getSon(1));

            // Create a list of partial likelihood update operations.
            // The order is [dest, destScaling, sourceScaling, source1, matrix1, source2, matrix2].
            // Possible TODO: no scaling supported here. Should there be?
            operations.push_back(BeagleOperation(
                                 {buffer,         // Destination buffer
                                  BEAGLE_OP_NONE, // (output) scaling buffer index
                                  BEAGLE_OP_NONE, // (input) scaling buffer index
                                  child1_buffer,  // Index of first child partials buffer
                                  child1_buffer,  // Index of first child transition matrix
                                  child2_buffer,  // Index of second child partials buffer
                                  child2_buffer})); // Index of second child transition matrix
            node_indices.push_back(child1_buffer);
            branch_lengths.push_back(n->getSon(0)->getDistanceToFather());
            node_indices.push_back(child2_buffer);
            branch_lengths.push_back(n->getSon(1)->getDistanceToFather());
        }
        distalNodeState[n] = hash_node(n);
    }

    updateTransitionsPartials(operations, branch_lengths, node_indices, nBuffers_);
    accumulateScaleFactors(operations, nBuffers_);
}

void BeagleTreeLikelihood::calculateProximalPartials()
{
    verifyInitialized();

    const std::vector<const bpp::Node*> preorder_nodes = preorder_find_changed(*tree, proxNodeState);

    // For tracking BEAGLE operations
    std::vector<BeagleOperation> operations;
    std::vector<int> node_indices;      // probability indices
    std::vector<double> branch_lengths;

    // Special handling for the two descendants of the root
    assert(tree->getRootNode()->getNumberOfSons() == 2);
    for(size_t i = 0; i < tree->getRootNode()->getNumberOfSons(); i++) {
        const bpp::Node* node = tree->getRootNode()->getSon(i);
        const bpp::Node* sibling = siblings(node)[0];
        proxNodeBuffer[node] = distalNodeBuffer[sibling];
        proxNodeState[node] = hash_node(node);
    }

    // Traverse internal nodes in preorder, adding BeagleOperations to update each
    for(const bpp::Node* n : preorder_nodes) {
        if(n->isLeaf() || n == tree->getRootNode())
            continue;
        assert(n->getNumberOfSons() == 2);
        // The distal likelihood for this node should already be calculated.
        assert(distalNodeBuffer.find(n) != distalNodeBuffer.end());

        const int parent_buffer = proxNodeBuffer.at(n);
        for(size_t i = 0; i < 2; ++i) {
            const bpp::Node* son = n->getSon(i);
            const bpp::Node* sibling = siblings(son).at(0);

            assert(distalNodeBuffer.find(sibling) != distalNodeBuffer.end());
            assert(proxNodeBuffer.find(son) != proxNodeBuffer.end());
            const int buffer = proxNodeBuffer.at(son);

            int sibling_buffer = distalNodeBuffer.at(sibling);

            operations.push_back(BeagleOperation(
                                 {buffer,            // Destination buffer
                                  BEAGLE_OP_NONE,    // (output) scaling buffer index
                                  BEAGLE_OP_NONE,    // (input) scaling buffer index
                                  parent_buffer,     // Index of first child partials buffer
                                  parent_buffer,     // Index of first child transition matrix
                                  sibling_buffer,    // Index of second child partials buffer
                                  sibling_buffer})); // Index of second child transition matrix
            node_indices.push_back(parent_buffer);
            branch_lengths.push_back(son->getDistanceToFather());
            node_indices.push_back(sibling_buffer);
            branch_lengths.push_back(sibling->getDistanceToFather());
            proxNodeState[son] = hash_node(son);
        }
    }

    updateTransitionsPartials(operations, branch_lengths, node_indices, nBuffers_ + 1);

    // No scale factor accumulations: we'll only use the proximal values for guided proposals
    // accumulate_scale_factors(operations, n_buffers + 1);
}

void BeagleTreeLikelihood::updateTransitionsPartials(const std::vector<BeagleOperation>& operations,
                                                         const std::vector<double>& branch_lengths,
                                                         const std::vector<int>& node_indices,
                                                         const int scaling_buffer)
{
    assert(branch_lengths.size() == node_indices.size());
    assert(branch_lengths.size() == 2 * operations.size());

    // Register topology, branch lengths; update transition matrices.
    beagle_check(beagleUpdateTransitionMatrices(beagleInstance_,        // instance
                                                0,                      // eigenIndex
                                                node_indices.data(),    // probabilityIndices
                                                NULL,                   // firstDerivativeIndices
                                                NULL,                   // secondDerivativeIndices
                                                branch_lengths.data(),  // edgeLengths
                                                node_indices.size()));   // count

    // Update partials for all traversed nodes
    beagle_check(beagleUpdatePartials(beagleInstance_, operations.data(), operations.size(), scaling_buffer));
}

std::vector<BeagleTreeLikelihood::NodePartials> BeagleTreeLikelihood::get_mid_edge_partials()
{
    calculateDistalPartials();
    calculateProximalPartials();

    std::vector<BeagleTreeLikelihood::NodePartials> result;

    const int buffer = scratchBuffers_[0];

    // Only calculate one mid-edge partial for the edge containing the root
    for(bpp::Node* node : online_available_edges(*tree))
    {
        assert(proxNodeBuffer.count(node) > 0);
        int prox_buffer = proxNodeBuffer.at(node);
        int dist_buffer = distalNodeBuffer.at(node);
        double d = node->getDistanceToFather();
        // Special handling for root node - add distance to other child
        if(node->getFather() == tree->getRootNode())
            d += siblings(node)[0]->getDistanceToFather();
        const double mid = d / 2;

        // Current partials should be set up
        std::vector<double> partials(partialLength());

        const std::vector<BeagleOperation> operations{
            BeagleOperation({buffer,         // Destination buffer
                            BEAGLE_OP_NONE,  // (output) scaling buffer index
                            BEAGLE_OP_NONE,  // (input) scaling buffer index
                            prox_buffer,     // Index of first child partials buffer
                            prox_buffer,     // Index of first child transition matrix
                            dist_buffer,     // Index of second child partials buffer
                            dist_buffer})};  // Index of second child transition matrix

        const std::vector<double> branch_lengths{mid,mid};
        const std::vector<int> node_indices{prox_buffer,dist_buffer};

        updateTransitionsPartials(operations, branch_lengths, node_indices, buffer);
        accumulateScaleFactors(operations, buffer);

        result.emplace_back(node, LikelihoodVector(nRates_, nSites_, nStates_));
        // TODO: Scale here?
        beagle_check(beagleGetPartials(beagleInstance_, buffer, BEAGLE_OP_NONE, result.back().second.get().data()));
    }
    return result;
}

LikelihoodVector BeagleTreeLikelihood::getDistalPartials(const bpp::Node* node)
{
    calculateDistalPartials();
    LikelihoodVector result(nRates_, nSites_, nStates_);
    const int buffer = distalNodeBuffer.at(node);
    beagle_check(beagleGetPartials(beagleInstance_, buffer, BEAGLE_OP_NONE, result.data()));
    return result;
}

LikelihoodVector BeagleTreeLikelihood::getProximalPartials(const bpp::Node* node)
{
    calculateProximalPartials();
    LikelihoodVector result(nRates_, nSites_, nStates_);
    const int buffer = proxNodeBuffer.at(node);
    beagle_check(beagleGetPartials(beagleInstance_, buffer, BEAGLE_OP_NONE, result.data()));
    return result;
}

LikelihoodVector BeagleTreeLikelihood::getLeafPartials(const std::string& name)
{
    const int buffer = leafBuffer.at(name);
    LikelihoodVector result(nRates_, nSites_, nStates_);
    beagle_check(beagleGetPartials(beagleInstance_, buffer, BEAGLE_OP_NONE, result.data()));
    return result;
}

int BeagleTreeLikelihood::getDistalBuffer(const bpp::Node* node)
{
    return distalNodeBuffer.at(node);
}

int BeagleTreeLikelihood::getProximalBuffer(const bpp::Node* node)
{
    return proxNodeBuffer.at(node);
}

int BeagleTreeLikelihood::getLeafBuffer(const std::string& name)
{
    return leafBuffer.at(name);
}

void BeagleTreeLikelihood::invalidate(const bpp::Node* node)
{
    distalNodeState.erase(node);
    proxNodeState.erase(node);
}

void BeagleTreeLikelihood::invalidateAll()
{
    distalNodeState.clear();
    proxNodeState.clear();
}

void BeagleTreeLikelihood::accumulateScaleFactors(const std::vector<BeagleOperation>& operations,
                                                      const int scale_buffer)
{
    // Calculate marginal log-likelihood scaling for each node
    const std::unique_ptr<int[]> scale_indices(new int[operations.size()]);
    for(size_t i = 0; i < operations.size(); i++) {
        scale_indices[i] = operations[i].destinationPartials;
    }
    beagle_check(beagleAccumulateScaleFactors(beagleInstance_, scale_indices.get(), operations.size(), scale_buffer));
}

double BeagleTreeLikelihood::calculateLogLikelihood()
{
    calculateDistalPartials();
    int root_buffer = distalNodeBuffer.at(tree->getRootNode());

    // Calculate root log likelihood
    const int category_weight_index = 0;
    const int state_frequency_index = 0;
    const int scaling_indices = nBuffers_;
    double log_likelihood;
    beagle_check(beagleCalculateRootLogLikelihoods(beagleInstance_,
                                                   &root_buffer,
                                                   &category_weight_index,
                                                   &state_frequency_index,
                                                   &scaling_indices,
                                                   1,
                                                   &log_likelihood));
    return log_likelihood;
}

void BeagleTreeLikelihood::verifyInitialized() const
{
    //if(tree == nullptr)
        //throw std::runtime_error("NULL tree");
    if(beagleInstance_ < 0)
        throw std::runtime_error("BEAGLE instance not initialized.");
}

}} // Namespace

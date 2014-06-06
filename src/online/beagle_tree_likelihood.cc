/// \file beagle_tree_likelihood.cc

#include "beagle_tree_likelihood.h"
#include "bpp_shim.h"
#include "likelihood_vector.h"
#include "util.h"
#include "online_util.h"

#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Seq/Container/SiteContainer.h>

#include <boost/graph/depth_first_search.hpp>

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

/// Mix-in for a depth first visitor that only visits nodes reachable from a given vertex root.
template<typename TGraph>
class SingleComponentMixIn : public boost::default_dfs_visitor
{
public:
    using TVertex = typename boost::graph_traits<TGraph>::vertex_descriptor;
    using TEdge = typename boost::graph_traits<TGraph>::edge_descriptor;
    using TEdgeIterator = typename boost::graph_traits<TGraph>::out_edge_iterator;

    explicit SingleComponentMixIn(const TVertex root)
    {
        inComponent[root] = true;
    }

    void examine_edge(TEdge edge, TGraph graph)
    {
        assert(inComponent.count(boost::target(edge, graph)) == 0 &&
               "Traversed an edge twice");
        inComponent[boost::target(edge, graph)] = inComponent[boost::source(edge, graph)];
    }

    bool operator()(TVertex vertex, TGraph)
    {
        return inComponent[vertex];
    }

    bool in_component(const TVertex vertex) const { return inComponent[vertex]; };
private:
    typename std::unordered_map<TVertex, bool> inComponent;
};

/// \brief Updates the hashes at all nodes reachable from the root, marks nodes and predecessors dirty when hash changes.
///
/// Run before BeagleUpdatePartialsVisitor
template<typename TGraph>
class BeagleMarkDirtyVisitor : public SingleComponentMixIn<TGraph>
{
public:
    using TVertex = typename boost::graph_traits<TGraph>::vertex_descriptor;
    using TEdge = typename boost::graph_traits<TGraph>::edge_descriptor;
    using TEdgeIterator = typename boost::graph_traits<TGraph>::out_edge_iterator;

    using SingleComponentMixIn<TGraph>::SingleComponentMixIn;

    void finish_edge(TEdge edge, TGraph graph)
    {
        const TVertex vertex = boost::target(edge, graph);
        if(in_component(vertex)) {
            // Rehash the node
            const double d = graph[edge];
            std::hash<TVertex> h;

            // Vertex ID
            size_t hash = h(vertex);
            // Edge length
            hash_combine(hash, d);
            // Child IDs
            TEdgeIterator it, end;
            std::tie(it, end) = boost::out_edges(vertex, graph);
            for(; it != end; it++) {
                if(graph[boost::target(*it, graph)].dirty) {
                    graph[vertex].dirty = true;
                }
                hash_combine(hash, boost::target(*it, graph));
            }
            if(hash != graph[vertex].hash || graph[vertex].dirty) {
                graph[vertex].hash = hash;
                graph[vertex].dirty = true;
            }
        }
    }
};

/// \brief Visit buffers in postorder, building lists of operations that need to be performed to bring the tree
/// partials up to date.
///
/// Example usage:
///
///     boost::vector_property_map<boost::default_color_type> colorVec(boost::num_vertices(graph));
///     boost::depth_first_visit(graph, vertex, visitor, colorVec, visitor);
///
/// This class should be run *after* BeagleMarkDirtyVisitor
template<typename TGraph>
class BeagleUpdatePartialsVisitor : public SingleComponentMixIn<TGraph>
{
public:
    using TVertex = typename boost::graph_traits<TGraph>::vertex_descriptor;
    using TEdge = typename boost::graph_traits<TGraph>::edge_descriptor;
    using TEdgeIterator = typename boost::graph_traits<TGraph>::out_edge_iterator;

    using SingleComponentMixIn<TGraph>::SingleComponentMixIn;

    void finish_vertex(TVertex vertex, TGraph graph)
    {
        TEdgeIterator it, end;
        std::tie(it, end) = boost::out_edges(vertex, graph);
        if(it == end) {
            // leaf
            return;
        }

        if(graph[vertex].dirty) {
            std::vector<int> targets;
            for(; it != end; ++it) {
                double dist = graph[*it];
                double target = graph[boost::target(*it, graph)].buffer;
                targets.push_back(target);
                branchLengths.push_back(dist);
            }

            assert(targets.size() == 2 && "Unexpected target size");
            const int buffer = graph[vertex].buffer;
            operations.push_back(BeagleOperation(
                                 {buffer,            // Destination buffer
                                  BEAGLE_OP_NONE,    // (output) scaling buffer index
                                  BEAGLE_OP_NONE,    // (input) scaling buffer index
                                  targets[0],        // Index of first child partials buffer
                                  targets[0],        // Index of first child transition matrix
                                  targets[1],        // Index of second child partials buffer
                                  targets[1]}));     // Index of second child transition matrix
            nodeIndices.insert(nodeIndices.end(), targets.begin(), targets.end());

            // Update hash
           graph[vertex].dirty = false;
        }
    }

    std::vector<BeagleOperation> operations;
    std::vector<int> nodeIndices;
    std::vector<double> branchLengths;

};

template<typename TGraph>
struct BeagleScaleFactorVisitor : public SingleComponentMixIn<TGraph>
{
    using TVertex = typename boost::graph_traits<TGraph>::vertex_descriptor;
    using TEdge = typename boost::graph_traits<TGraph>::edge_descriptor;

    using SingleComponentMixIn<TGraph>::SingleComponentMixIn;

    void finish_vertex(TVertex vertex, TGraph graph)
    {
        using TEdgeIterator = typename boost::graph_traits<TGraph>::out_edge_iterator;
        TEdgeIterator it, end;
        std::tie(it, end) = boost::out_edges(vertex, graph);
        // Leaf nodes have no out edges, and do not need to be scaled.
        if(it != end) {
            buffers.push_back(graph[vertex].buffer);
        }
    }

    std::unordered_map<TVertex, bool> inComponent;
    std::vector<int> buffers;
};


BeagleTreeLikelihood::BeagleTreeLikelihood(const bpp::SiteContainer& sites,
                                           const bpp::SubstitutionModel& model,
                                           const bpp::DiscreteDistribution& rateDist,
                                           const size_t nScratchBuffers) :
    beagleInstance_(-1),
    nSites_(sites.getNumberOfSites()),
    nStates_(model.getNumberOfStates()),
    nRates_(rateDist.getNumberOfCategories()),
    nSeqs_(sites.getNumberOfSequences()),
    // Allocate three buffers for each node in the tree (to store distal, proximal vectors, mid-edge vectors)
    // plus `scratch_buffer_count` BONUS buffers
    nBuffers_((2 * nSeqs_ - 1) * 3 + nScratchBuffers),
    rateDist(&rateDist),
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

    // Fill available buffers
    for(size_t i = 0; i < nBuffers_; i++)
        availableBuffers.push(nBuffers_ - 1 - i);

    // Load tips
    for(size_t i = 0; i < nSeqs_; i++)
        registerLeaf(sites.getSequence(i));

    // Weight all sites equally -
    // for online inference, we don't want to compress sites.
    std::vector<double> pattern_weights(sites.getNumberOfSites(), 1.0);
    beagle_check(beagleSetPatternWeights(beagleInstance_, pattern_weights.data()));
}

BeagleTreeLikelihood::~BeagleTreeLikelihood()
{
    if(beagleInstance_ >= 0)
        beagleFinalizeInstance(beagleInstance_);
}

int BeagleTreeLikelihood::getFreeBuffer()
{
    assert(!availableBuffers.empty());
    const int buffer = availableBuffers.top();
    availableBuffers.pop();

    assert(bufferMap.find(buffer) == bufferMap.end() && "used buffer in available buffers.");

    addBufferToGraph(BeagleTreeLikelihood::VertexInfo{buffer, 0, true});

    return buffer;
}

BeagleBuffer BeagleTreeLikelihood::borrowBuffer()
{
    assert(!availableBuffers.empty());
    return BeagleBuffer(this);
}

void BeagleTreeLikelihood::returnBuffer(const int buffer, const bool check)
{
    assert(buffer < nBuffers_);
    auto it = bufferMap.find(buffer);
    if(check) {
        assert(it != bufferMap.end() && "Tried to return unknown buffer!");
    }
    if(it != bufferMap.end()) {
        boost::remove_vertex(it->second, graph);
        bufferMap.erase(it);
        availableBuffers.push(buffer);
    }
}

void BeagleTreeLikelihood::initialize(const bpp::SubstitutionModel& model,
                                      const bpp::DiscreteDistribution& rateDist,
                                      bpp::TreeTemplate<bpp::Node>& tree)
{
    this->rateDist = &rateDist;
    this->model = &model;
    this->tree = &tree;
    verifyInitialized();

    // Clear buffer maps
    distalNodeBuffer.clear();
    proxNodeBuffer.clear();
    midEdgeNodeBuffer.clear();
    graph = TGraph();
    invalidateAll();

    // Return all buffers that aren't associated with a leaf.
    std::vector<bool> isNonLeafBuffer(nBuffers_, true);
    for(const auto& p : leafBuffer) {
        isNonLeafBuffer[graph[p.second].buffer] = false;
    }
    assert(std::count(isNonLeafBuffer.begin(), isNonLeafBuffer.end(), false) == leafBuffer.size());
    for(size_t i = 0; i < nBuffers_; i++) {
        if(isNonLeafBuffer[i]) {
            returnBuffer(i, false);
        }
    }

    loadRateDistribution(rateDist);
    loadSubstitutionModel(model);

    // Slide root to far right side of root branch
    tree.getRootNode()->getSon(0)->setDistanceToFather(tree.getRootNode()->getSon(0)->getDistanceToFather() +
                                                       tree.getRootNode()->getSon(1)->getDistanceToFather());
    tree.getRootNode()->getSon(1)->setDistanceToFather(0.0);

    // Fill buffer maps

    // Mid-edge buffer
    allocateDistalBuffers();
    allocateProximalBuffers();
    allocateMidEdgeBuffers();
    buildBufferDependencyGraph();
}

BeagleTreeLikelihood::TVertex BeagleTreeLikelihood::addBufferToGraph(const VertexInfo& info)
{
    if(bufferMap.find(info.buffer) != bufferMap.end())
        throw std::runtime_error("Buffer " + std::to_string(info.buffer) + " is already in graph.");
    TVertex vertex = boost::add_vertex(info, graph);
    bufferMap[info.buffer] = vertex;
    return vertex;
}

void BeagleTreeLikelihood::allocateDistalBuffers()
{
    const std::vector<bpp::Node*> nodes = tree->getNodes();
    for(const bpp::Node* n : nodes) {
        if(n->isLeaf()) {
            const std::string& name = n->getName();
            assert(leafBuffer.count(name) > 0);
            distalNodeBuffer[n] = leafBuffer.at(name);
        } else {
            assert(n->getNumberOfSons() == 2);
            assert(distalNodeBuffer.find(n) == distalNodeBuffer.end());
            distalNodeBuffer[n] = bufferMap.at(getFreeBuffer());
        }
    }
}

void BeagleTreeLikelihood::allocateProximalBuffers()
{
    // Special handling at the root - proximal buffers for each child is the *distal* buffer for its sibling.
    const bpp::Node* root = tree->getRootNode();
    for(size_t i = 0; i < root->getNumberOfSons(); i++) {
        const bpp::Node* n = root->getSon(i);
        proxNodeBuffer[n] = distalNodeBuffer.at(siblings(n)[0]);
    }
    for(const bpp::Node* n : tree->getNodes()) {
        if(n->isLeaf() || n == root)
            continue;
        assert(n->getNumberOfSons() == 2);
        for(size_t i = 0; i < 2; ++i) {
            const bpp::Node* son = n->getSon(i);
            assert(distalNodeBuffer.find(son) != distalNodeBuffer.end());
            assert(proxNodeBuffer.find(son) == proxNodeBuffer.end());
            proxNodeBuffer[son] = bufferMap.at(getFreeBuffer());
        }
    }
}

void BeagleTreeLikelihood::allocateMidEdgeBuffers()
{
    for(const bpp::Node* n : onlineAvailableEdges(*tree)) {
        assert(n != nullptr);
        midEdgeNodeBuffer[n] = bufferMap.at(getFreeBuffer());
    }
}

void BeagleTreeLikelihood::buildBufferDependencyGraph()
{
    // Distal
    // This is easy - each node just depends on the distal buffers of its children.
    for(const bpp::Node* n : tree->getNodes()) {
        if(n->isLeaf())
            continue;
        assert(n->getNumberOfSons() == 2);
        for(size_t i = 0; i < n->getNumberOfSons(); i++) {
            const bpp::Node* son = n->getSon(i);
            boost::add_edge(distalNodeBuffer.at(n),
                            distalNodeBuffer.at(son),
                            son->getDistanceToFather(),
                            graph);
        }
    }

    // Proximal
    // Here, we traverse nodes, adding proximal nodes for the sons of each non-leaf, non-root node.
    for(const bpp::Node* parent : tree->getNodes()) {
        if(parent->isLeaf() || parent == tree->getRootNode())
            continue;
        assert(parent->getNumberOfSons() == 2);
        // The distal buffer for this node should already be calculated.
        assert(distalNodeBuffer.find(parent) != distalNodeBuffer.end());

        const TVertex parentVertex = proxNodeBuffer.at(parent);

        for(size_t i = 0; i < parent->getNumberOfSons(); ++i) {
            const bpp::Node* son = parent->getSon(i);
            const bpp::Node* sibling = siblings(son).at(0);

            assert(distalNodeBuffer.find(sibling) != distalNodeBuffer.end());
            assert(proxNodeBuffer.find(son) != proxNodeBuffer.end());
            const TVertex vertex = proxNodeBuffer.at(son);

            const TVertex siblingVertex = distalNodeBuffer.at(sibling);

            boost::add_edge(vertex, siblingVertex, sibling->getDistanceToFather(), graph);
            double parentDist = parent->getDistanceToFather();
            if(parent->getFather() == tree->getRootNode())
                parentDist += siblings(parent)[0]->getDistanceToFather();
            boost::add_edge(vertex, parentVertex, parentDist, graph);
        }
    }

    // Mid-Edge buffers - depend on the proximal and distal buffers of the edge.
    for(const bpp::Node* n : onlineAvailableEdges(*tree)) {
        const TVertex prox = proxNodeBuffer.at(n),
                      distal = distalNodeBuffer.at(n),
                      midEdge = midEdgeNodeBuffer.at(n);
        double d = n->getDistanceToFather();
        // Special handling for root node - distance should be sum of branches below root
        if(n->getFather() == tree->getRootNode())
            d += siblings(n)[0]->getDistanceToFather();
        const double mid = d / 2;
        boost::add_edge(midEdge, prox, mid, graph);
        boost::add_edge(midEdge, distal, mid, graph);
    }
}

size_t BeagleTreeLikelihood::registerLeaf(const bpp::Sequence& sequence)
{
    verifyInitialized();
    const int buffer = getFreeBuffer();
    graph[bufferMap.at(buffer)].dirty = false;
    if(leafBuffer.count(sequence.getName()) > 0)
        throw std::runtime_error("Duplicate sequence name: " + sequence.getName());

    leafBuffer[sequence.getName()] = bufferMap.at(buffer);
    const std::vector<double> seq_partials = get_partials(sequence, *model, nRates_);
    assert(seq_partials.size() == sequence.size() * nStates_ * nRates_);
    beagle_check(beagleSetPartials(beagleInstance_, buffer, seq_partials.data()));
    return buffer;
}

void BeagleTreeLikelihood::loadSubstitutionModel(const bpp::SubstitutionModel& model)
{
    verifyInitialized();
    // Eigendecomposition
    std::vector<double> evec(nStates_ * nStates_),
                        ivec(nStates_ * nStates_),
                        eval(nStates_);
    blit_matrix_to_array(ivec.data(), model.getRowLeftEigenVectors());     // inverse eigenvectors
    blit_matrix_to_array(evec.data(), model.getColumnRightEigenVectors()); // eigenvectors
    blit_vector_to_array(eval.data(), model.getEigenValues());
    beagle_check(beagleSetEigenDecomposition(beagleInstance_, 0, evec.data(), ivec.data(), eval.data()));
    // State frequencies
    beagle_check(beagleSetStateFrequencies(beagleInstance_, 0, model.getFrequencies().data()));
}

void BeagleTreeLikelihood::loadRateDistribution(const bpp::DiscreteDistribution& rateDist)
{
    assert(rateDist.getNumberOfCategories() == nRates_ &&
           "Unexpected rate category count");
    verifyInitialized();
    const std::vector<double>& categories = rateDist.getCategories();
    const std::vector<double>& weights = rateDist.getProbabilities();

    beagle_check(beagleSetCategoryRates(beagleInstance_, categories.data()));
    beagle_check(beagleSetCategoryWeights(beagleInstance_, 0, weights.data()));
}

void BeagleTreeLikelihood::updateTransitionsPartials(const std::vector<BeagleOperation>& operations,
                                                     const std::vector<double>& branchLengths,
                                                     const std::vector<int>& nodeIndices,
                                                     const int scalingBuffer)
{
    assert(branchLengths.size() == nodeIndices.size());
    assert(branchLengths.size() == 2 * operations.size());

    // Register topology, branch lengths; update transition matrices.
    beagle_check(beagleUpdateTransitionMatrices(beagleInstance_,        // instance
                                                0,                      // eigenIndex
                                                nodeIndices.data(),     // probabilityIndices
                                                NULL,                   // firstDerivativeIndices
                                                NULL,                   // secondDerivativeIndices
                                                branchLengths.data(),   // edgeLengths
                                                nodeIndices.size()));   // count

    // Update partials for all traversed nodes
    beagle_check(beagleUpdatePartials(beagleInstance_, operations.data(),
                                      operations.size(), scalingBuffer));
}

void BeagleTreeLikelihood::updateTransitionsPartials(const TVertex vertex)
{
    // First, update dirty/clean status
    {
        BeagleMarkDirtyVisitor<TGraph> visitor(vertex);
        boost::vector_property_map<boost::default_color_type> colorVec(boost::num_vertices(graph));
        boost::depth_first_visit(graph, vertex, visitor, colorVec, visitor);
    }
    // Now update the partials of any dirty nodes
    BeagleUpdatePartialsVisitor<TGraph> visitor(vertex);
    boost::vector_property_map<boost::default_color_type> colorVec(boost::num_vertices(graph));
    boost::depth_first_visit(graph, vertex, visitor, colorVec, visitor);
    updateTransitionsPartials(visitor.operations, visitor.branchLengths, visitor.nodeIndices, BEAGLE_OP_NONE);
}

std::vector<BeagleTreeLikelihood::NodePartials> BeagleTreeLikelihood::getMidEdgePartials()
{
    std::vector<bpp::Node*> nodes = onlineAvailableEdges(*tree);
    std::vector<BeagleTreeLikelihood::NodePartials> result;
    result.reserve(nodes.size());

    // Only calculate one mid-edge partial for the edge containing the root
    // onlineAvailableEdges skips edge to the right of the root
    for(bpp::Node* node : nodes)
    {
        TVertex midEdgeVertex = midEdgeNodeBuffer.at(node);
        updateTransitionsPartials(midEdgeVertex);
        assert(!graph[midEdgeVertex].dirty && "Vertex should be clean");
        result.emplace_back(node, graph[midEdgeVertex].buffer);
    }
    return result;
}

std::vector<double> BeagleTreeLikelihood::calculateAttachmentLikelihood(const std::string& leafName,
                                                                        const bpp::Node* node,
                                                                        const double distalLength,
                                                                        const std::vector<double>& pendantBranchLengths)
{
    if(!leafBuffer.count(leafName))
        throw std::runtime_error("Unknown leaf: " + leafName);

    const int leafBuf = leafBuffer[leafName];

    const BeagleBuffer b = borrowBuffer();
    const TVertex vert = bufferMap.at(b.value());

    // TODO: CHECK
    double edgeLength = node->getDistanceToFather();
    if(node->getFather() == tree->getRootNode())
        edgeLength += siblings(node)[0]->getDistanceToFather();

    if(distalLength > edgeLength)
        throw std::runtime_error("Invalid distal length!");

    const TVertex dist = distalNodeBuffer.at(node),
                  prox = proxNodeBuffer.at(node);
    addDependencies(vert, prox, distalLength, dist, edgeLength - distalLength);
    updateTransitionsPartials(vert);

    std::vector<double> result;
    result.reserve(pendantBranchLengths.size());
    for(const double pendant : pendantBranchLengths)
        result.push_back(logDot(b.value(), leafBuf, pendant));
    return result;
}


std::vector<std::vector<double>> BeagleTreeLikelihood::calculateAttachmentLikelihoods(const std::string& leafName,
                                                                                      const std::vector<BeagleTreeLikelihood::AttachmentLocation>& attachmentLocations,
                                                                                      const std::vector<double> pendantBranchLengths)
{
    if(!leafBuffer.count(leafName))
        throw std::runtime_error("Unknown leaf: " + leafName);

    std::vector<std::vector<double>> result;
    result.reserve(attachmentLocations.size());

    for(auto& loc : attachmentLocations) {
        result.push_back(calculateAttachmentLikelihood(leafName, loc.first, loc.second, pendantBranchLengths));
    }

    return result;
}

LikelihoodVector BeagleTreeLikelihood::getDistalPartials(const bpp::Node* node)
{
    LikelihoodVector result(nRates_, nSites_, nStates_);
    const TVertex vertex = distalNodeBuffer.at(node);
    updateTransitionsPartials(vertex);
    beagle_check(beagleGetPartials(beagleInstance_, graph[vertex].buffer, BEAGLE_OP_NONE, result.data()));
    return result;
}

LikelihoodVector BeagleTreeLikelihood::getProximalPartials(const bpp::Node* node)
{
    LikelihoodVector result(nRates_, nSites_, nStates_);
    const TVertex vertex = proxNodeBuffer.at(node);
    updateTransitionsPartials(vertex);
    beagle_check(beagleGetPartials(beagleInstance_, graph[vertex].buffer, BEAGLE_OP_NONE, result.data()));
    return result;
}

LikelihoodVector BeagleTreeLikelihood::getLeafPartials(const std::string& name) const
{
    const int buffer = graph[leafBuffer.at(name)].buffer;
    LikelihoodVector result(nRates_, nSites_, nStates_);
    beagle_check(beagleGetPartials(beagleInstance_, buffer, BEAGLE_OP_NONE, result.data()));
    return result;
}

int BeagleTreeLikelihood::getDistalBuffer(const bpp::Node* node) const
{
    return graph[distalNodeBuffer.at(node)].buffer;
}

int BeagleTreeLikelihood::getProximalBuffer(const bpp::Node* node) const
{
    return graph[proxNodeBuffer.at(node)].buffer;
}

int BeagleTreeLikelihood::getMidEdgeBuffer(const bpp::Node* node) const
{
    return graph[midEdgeNodeBuffer.at(node)].buffer;
}

int BeagleTreeLikelihood::getLeafBuffer(const std::string& name) const
{
    return graph[leafBuffer.at(name)].buffer;
}

void BeagleTreeLikelihood::invalidate(const bpp::Node* node)
{
    graph[distalNodeBuffer.at(node)].dirty = true;
    graph[distalNodeBuffer.at(node)].hash = 0;
    graph[proxNodeBuffer.at(node)].dirty = true;
    graph[proxNodeBuffer.at(node)].hash = 0;
}

void BeagleTreeLikelihood::invalidateAll()
{
    using TIterator = boost::graph_traits<TGraph>::vertex_iterator;
    TIterator it, end;
    for(; it != end; ++it) {
        graph[*it].dirty = true;
        graph[*it].hash = 0;
    }
}

void BeagleTreeLikelihood::accumulateScaleFactors(const std::vector<BeagleOperation>& operations,
                                                  const int scaleBuffer)
{
    // Calculate marginal log-likelihood scaling for each node
    std::vector<int> scaleIndices(operations.size(), 0);
    for(size_t i = 0; i < operations.size(); i++) {
        scaleIndices[i] = operations[i].destinationPartials;
    }
    beagle_check(beagleAccumulateScaleFactors(beagleInstance_, scaleIndices.data(), operations.size(), scaleBuffer));
}

double BeagleTreeLikelihood::calculateLogLikelihood()
{
    TVertex root = distalNodeBuffer.at(tree->getRootNode());
    updateTransitionsPartials(root);
    return logLikelihood(graph[root].buffer);
}

void BeagleTreeLikelihood::verifyInitialized() const
{
    //if(tree == nullptr)
        //throw std::runtime_error("NULL tree");
    if(beagleInstance_ < 0)
        throw std::runtime_error("BEAGLE instance not initialized.");
}

void BeagleTreeLikelihood::toDot(std::ostream& out) const
{
    out << "digraph {\n";
    for(const bpp::Node* n : postorder(tree->getRootNode())) {
        const int nodeId = n->getId();
        std::string label = 'n' + std::to_string(nodeId);
        if(n->isLeaf())
            label += " [" + n->getName() + ']';
        else if (!n->hasFather())
            label += " [root]";
        out << nodeId << "[label=\"" << label << "\"];\n";
        if(n->hasFather())
            out << n->getFatherId() << " -> " << nodeId
                << "[label=" << n->getDistanceToFather() << "];\n";
    }

    const std::vector<bpp::Node*> nodes = tree->getNodes();
    // Distal buffer
    for(const bpp::Node* n : nodes) {
        const int distalBuffer = distalNodeBuffer.at(n);
        out << "b" << distalBuffer << "[shape=none];\n";
        out << "b" << distalBuffer << " -> " << n->getId() << "[color=blue];\n";
    }
    const bpp::Node* root = tree->getRootNode();
    for(size_t i = 0; i < root->getNumberOfSons(); i++) {
        const bpp::Node* n = root->getSon(i);
        const int proxBuffer = proxNodeBuffer.at(n);
        out << "b" << proxBuffer << "[shape=none];\n";
        out << "b" << proxBuffer << " -> " << n->getId() << "[color=red,style=dashed];\n";
    }
    for(const bpp::Node* n : nodes) {
        if(n->isLeaf() || n == root)
            continue;
        for(size_t i = 0; i < 2; ++i) {
            const bpp::Node* son = n->getSon(i);
            const int proxBuffer = proxNodeBuffer.at(son);
            out << "b" << proxBuffer << "[shape=none];\n";
            out << "b" << proxBuffer << " -> " << son->getId() << "[color=red,style=dashed];\n";
        }
    }

    for(const auto &p : midEdgeNodeBuffer) {
        const int distBuffer = distalNodeBuffer.at(p.first);
        const int proxBuffer = proxNodeBuffer.at(p.first);
        out << "b" << p.second << "[shape=none];\n";
        out << "b" << p.second << " -> b" << distBuffer << "[color=green,style=dotted]\n";
        out << "b" << p.second << " -> b" << proxBuffer << "[color=green,style=dotted]\n";
    }

    out << "}\n";
}

double BeagleTreeLikelihood::logDot(const std::vector<double>& v1,
                                    const std::vector<double>& v2,
                                    const double d)
{
    assert(v2.size() == partialLength() && "unexpected partial length");
    assert(freeBufferCount() >= 3);
    const BeagleBuffer b = borrowBuffer();
    beagleSetPartials(beagleInstance_, b.value(), v2.data());
    return logDot(v1, b.value(), d);
}

double BeagleTreeLikelihood::logDot(const std::vector<double>& v, const int buffer, const double d)
{
    assert(v.size() == partialLength() && "unexpected partial length");
    assert(freeBufferCount() >= 2);
    const BeagleBuffer b = borrowBuffer();
    beagleSetPartials(beagleInstance_, b.value(), v.data());

    return logDot(b.value(), buffer, d);
}

double BeagleTreeLikelihood::logDot(const int buffer1, const int buffer2, const double d)
{
    assert(buffer1 < nBuffers_ && buffer1 >= 0 && "Invalid buffer!");
    assert(buffer2 < nBuffers_ && buffer2 >= 0 && "Invalid buffer!");
    assert(freeBufferCount() >= 1);

    const BeagleBuffer b = borrowBuffer();
    const int scratchBuffer = b.value();
    assert(scratchBuffer != buffer1 && scratchBuffer != buffer2 &&
           "Reused buffer");

    addDependencies(bufferMap.at(scratchBuffer), bufferMap.at(buffer1), d, bufferMap.at(buffer2), 0);
    updateTransitionsPartials(bufferMap.at(scratchBuffer));

    return logLikelihood(scratchBuffer);
}

double BeagleTreeLikelihood::logLikelihood(const std::vector<double>& v)
{
    assert(v.size() == partialLength() && "unexpected partial length");
    const BeagleBuffer b = borrowBuffer();
    const int tmpBuffer = b.value();
    if(!(instanceDetails.flags & BEAGLE_FLAG_SCALING_AUTO))
        beagleResetScaleFactors(beagleInstance_, tmpBuffer);
    beagleSetPartials(beagleInstance_, tmpBuffer, v.data());
    return logLikelihood(tmpBuffer);
}

double BeagleTreeLikelihood::logLikelihood(const int buffer)
{
    const int category_weight_index = 0;
    const int state_frequency_index = 0;
    const int scalingIndex = BEAGLE_OP_NONE;
    double log_likelihood;

    // Update scale factors
    const TVertex vertex = bufferMap.at(buffer);
    BeagleScaleFactorVisitor<TGraph> visitor(vertex);
    boost::vector_property_map<boost::default_color_type> colorVec(boost::num_vertices(graph));
    boost::depth_first_visit(graph, vertex, visitor, colorVec, visitor);

    beagle_check(beagleAccumulateScaleFactors(beagleInstance_,
                                              visitor.buffers.data(),
                                              visitor.buffers.size(),
                                              BEAGLE_OP_NONE));

    beagle_check(beagleCalculateRootLogLikelihoods(beagleInstance_,
                                                   &buffer,
                                                   &category_weight_index,
                                                   &state_frequency_index,
                                                   &scalingIndex,
                                                   1,  // op count
                                                   &log_likelihood));
    return log_likelihood;
}

}} // Namespace

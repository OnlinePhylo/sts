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
#include <iostream>
#include <stack>
#include <stdexcept>
#include <unordered_set>


using sts::likelihood::blit_vector_to_array;
using sts::likelihood::blit_matrix_to_array;
using sts::likelihood::get_partials;
using sts::util::beagle_check;


namespace sts { namespace online {

/// From boost: add a hash of `v` to `seed`
template<typename T> void hash_combine(size_t& seed, const T& v)
{
    std::hash<T> h;
    seed ^= h(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

/// Mix-in for a depth first visitor that only visits nodes reachable from a given vertex root.
/// \see [Boost DFSVisitor documentation](http://www.boost.org/doc/libs/1_55_0/libs/graph/doc/DFSVisitor.html)
/// \note Visitors are passed by value to `boost::depth_first_visit`.
/// Hold references to anything needing to persist.
///
/// \note The input graph is passed as a *const reference or by value* - it cannot be modified.
///
/// Process:
/// 1. SingleComponentMixIn::start_vertex colors the root node.
/// 2. SingleComponentMixIn::examine_edge extends coloring from parent to child vertices before the child vertex is
///    visited.
/// 3. SingleComponentMixIn::operator() terminates the search at the first uncolored vertex encountered -
///    this is the first vertex *outside* the tree rooted at the root passed to
///    SingleComponentMixIn::start_vertex.
///
/// \tparam TGraph Graph type (see BeagleTreeLikelihood::TGraph)
template<typename TGraph>
class SingleComponentMixIn : public boost::default_dfs_visitor
{
public:
    using TVertex = typename boost::graph_traits<TGraph>::vertex_descriptor;
    using TEdge = typename boost::graph_traits<TGraph>::edge_descriptor;
    using TEdgeIterator = typename boost::graph_traits<TGraph>::out_edge_iterator;

    /// \brief Initialize the search by coloring `root` as a member of the component of interest,
    virtual void start_vertex(TVertex root, const TGraph&)
    {
        assert(inComponent.size() == 0);
        inComponent[root] = true;
    }

    /// \brief Extend the coloring of `boost::source(edge)` to `boost::target(edge)`.
    virtual void examine_edge(TEdge edge, const TGraph& graph)
    {
        if(inComponent.count(boost::target(edge, graph))) {
            throw std::runtime_error("Re-visiting node: " + std::to_string(boost::target(edge, graph)));
        }
        inComponent[boost::target(edge, graph)] = inComponent[boost::source(edge, graph)];
    }

    /// \brief Tree filter
    ///
    /// This prevents the visitor from examining multiple connected components of the graph.
    /// \see [boost::depth_first_visit](http://www.boost.org/doc/libs/1_55_0/libs/graph/doc/depth_first_visit.html).
    bool operator()(const TVertex vertex, const TGraph&) const
    {
        auto it = inComponent.find(vertex);
        if(it == inComponent.end())
            return false;
        return it->second;
    }

    /// \brief Throw an error when there are multiple paths to the same node
    void forward_or_cross_edge(TEdge e, const TGraph& g)
    {
        throw std::runtime_error("Forward / cross edge between " + std::to_string(boost::source(e, g)) + "->" + std::to_string(boost::target(e, g)));
    }

    /// \brief Throw an error when a cycle is detected.
    void back_edge(TEdge e, const TGraph& g)
    {
        throw std::runtime_error("Back edge between " + std::to_string(boost::source(e, g)) + "->" + std::to_string(boost::target(e, g)));
    }

    bool in_component(const TVertex vertex) const { return inComponent.at(vertex); };
private:
    typename std::unordered_map<TVertex, bool> inComponent;
};

/// \brief Updates the hashes at all nodes reachable from the root, marks nodes and predecessors dirty when hash changes.
///
/// Run before BeagleUpdatePartialsVisitor.
/// \note The visitor gets a const reference to the graph.
///       After running, the graph state must be updated via BeagleMarkDirtyVisitor::update_graph.
/// \see [Boost DFSVisitor](http://www.boost.org/doc/libs/1_55_0/libs/graph/doc/DFSVisitor.html)
/// \tparam TGraph Graph type (see BeagleTreeLikelihood::TGraph)
template<typename TGraph>
class BeagleMarkDirtyVisitor : public SingleComponentMixIn<TGraph>
{
public:
    using TVertex = typename boost::graph_traits<TGraph>::vertex_descriptor;
    using TEdge = typename boost::graph_traits<TGraph>::edge_descriptor;
    using TEdgeIterator = typename boost::graph_traits<TGraph>::out_edge_iterator;
    using TVertexInfo = typename boost:: vertex_bundle_type<TGraph>::type;

    explicit BeagleMarkDirtyVisitor(std::unordered_map<TVertex, TVertexInfo>& result)
        : vertexInfo(result) {};

    virtual void start_vertex(TVertex root, const TGraph& graph)
    {
        SingleComponentMixIn<TGraph>::start_vertex(root, graph);
        vertexInfo.reserve(boost::num_vertices(graph));
    }

    /// \brief Rehash `vertex`
    ///
    /// When we finish a vertex, its target's children have already been visited.
    /// A node is then dirty when a) either of its children are dirty, or b) its state changed (by adding/removing a
    /// child, or edge lengths to the children changed.
    ///
    /// Leaves are always clean.
    void finish_vertex(const TVertex vertex, const TGraph& graph)
    {
        if(!this->in_component(vertex))
            return;

        vertexInfo[vertex] = graph[vertex];
        TVertexInfo& info = vertexInfo[vertex];
        if(graph[vertex].leaf) {
            assert(graph[vertex].dirty == false && "leaves should be clean.");
            assert(graph[vertex].hash == 0 && "leaves should never be hashed.");
            return;
        } else {
            // Rehash the node
            std::hash<TVertex> h;

            // Vertex ID
            size_t hash = h(vertex);

            // Children
            TEdgeIterator it, end;
            std::tie(it, end) = boost::out_edges(vertex, graph);
            for(; it != end; it++) {
                if(vertexInfo.at(boost::target(*it, graph)).dirty) {
                    info.dirty = true;
                }
                // Child ID
                hash_combine(hash, boost::target(*it, graph));
                // Distance to child
                hash_combine(hash, graph[*it]);
            }
            if(hash != info.hash || info.dirty) {
                info.hash = hash;
                info.dirty = true;
            }
        }
    }

    void update_graph(TGraph& graph) const
    {
        for(auto& p : vertexInfo) {
            graph[p.first] = p.second;
        }
    }
private:
    std::unordered_map<TVertex, TVertexInfo>& vertexInfo;
};

/// \brief Visit buffers in postorder, building lists of operations that need to be performed to bring the tree
/// partials up to date.
///
/// Example usage:
///
///     boost::vector_property_map<boost::default_color_type> colorVec(boost::num_vertices(graph));
///     boost::depth_first_visit(graph, vertex, visitor, colorVec, visitor);
///
/// This class should be run *after* BeagleMarkDirtyVisitor.
/// After running depth_first_search, mark nodes clean via BeagleUpdatePartialsVisitor::update_graph.
/// \see [Boost DFSVisitor](http://www.boost.org/doc/libs/1_55_0/libs/graph/doc/DFSVisitor.html)
/// \tparam TGraph Graph type (see BeagleTreeLikelihood::TGraph)
template<typename TGraph>
class BeagleUpdatePartialsVisitor : public SingleComponentMixIn<TGraph>
{
public:
    using TVertex = typename boost::graph_traits<TGraph>::vertex_descriptor;
    using TEdge = typename boost::graph_traits<TGraph>::edge_descriptor;
    using TEdgeIterator = typename boost::graph_traits<TGraph>::out_edge_iterator;
    using TVertexInfo = typename boost:: vertex_bundle_type<TGraph>::type;

    BeagleUpdatePartialsVisitor(std::vector<BeagleOperation>& operations,
                                std::vector<int>& nodeIndices,
                                std::vector<double>& branchLengths,
                                std::unordered_set<TVertex>& u) :
        operations(operations),
        nodeIndices(nodeIndices),
        branchLengths(branchLengths),
        updated(u) {}

    void finish_vertex(TVertex vertex, const TGraph& graph)
    {
        if(!this->in_component(vertex))
            return;

        TEdgeIterator it, end;
        std::tie(it, end) = boost::out_edges(vertex, graph);
        if(it == end) {
            // leaf
            return;
        }

        if(graph[vertex].dirty) {
            updated.insert(vertex);
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
        }
    }

    void update_graph(TGraph& graph) const
    {
        for(TVertex v : updated) {
            graph[v].dirty = false;
        }
    }
private:
    std::vector<BeagleOperation>& operations;
    std::vector<int>& nodeIndices;
    std::vector<double>& branchLengths;
    std::unordered_set<TVertex>& updated;
};

/// \brief Populates a vector scale buffers to accumulate via beagleAccumulateScaleFactors
template<typename TGraph>
class BeagleScaleFactorVisitor : public SingleComponentMixIn<TGraph>
{
public:
    using TVertex = typename boost::graph_traits<TGraph>::vertex_descriptor;
    using TEdge = typename boost::graph_traits<TGraph>::edge_descriptor;

    BeagleScaleFactorVisitor(std::vector<int>& buffers) :
        buffers(buffers) {};

    void finish_vertex(const TVertex vertex, const TGraph& graph)
    {
        if(!this->in_component(vertex))
            return;

        using TEdgeIterator = typename boost::graph_traits<TGraph>::out_edge_iterator;
        TEdgeIterator it, end;
        std::tie(it, end) = boost::out_edges(vertex, graph);
        // Leaf nodes have no out edges, and do not need to be scaled.
        if(it != end) {
            assert(!graph[vertex].leaf);
            buffers.push_back(graph[vertex].buffer);
        }
    }

private:
    std::vector<int>& buffers;
};


// Initialize static var
size_t BeagleTreeLikelihood::totalBeagleUpdateTransitionsCalls_ = 0;

BeagleTreeLikelihood::BeagleTreeLikelihood(const bpp::SiteContainer& sites,
                                           const bpp::SubstitutionModel& model,
                                           const bpp::DiscreteDistribution& rateDist,
                                           const size_t nScratchBuffers) :
    beagleInstance_(-1),
    nSites_(sites.getNumberOfSites()),
    nStates_(model.getNumberOfStates()),
    nRates_(rateDist.getNumberOfCategories()),
    nSeqs_(sites.getNumberOfSequences()),
    // Allocate two buffers for each node in the tree (to store distal, proximal vectors)
    // plus `scratch_buffer_count` BONUS buffers
    nBuffers_((2 * nSeqs_ - 1) * 2 + nScratchBuffers),
    nBeagleUpdateTransitionsCalls_(0),
    rateDist_(&rateDist),
    model_(&model),
    tree_(nullptr)
{
    assert(nRates_ >= 1);

    leafVertex_.reserve(nSeqs_);

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
            &instanceDetails_);
    if(beagleInstance_ < 0)
        beagle_check(beagleInstance_);

    // Fill available buffers
    for(size_t i = 0; i < nBuffers_; i++)
        availableBuffers_.push(nBuffers_ - 1 - i);

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
    assert(!availableBuffers_.empty());
    const int buffer = availableBuffers_.top();
    availableBuffers_.pop();

    assert(bufferMap_.find(buffer) == bufferMap_.end() && "used buffer in available buffers.");

    addBufferToGraph(BeagleTreeLikelihood::VertexInfo{buffer, 0, true, false});

    return buffer;
}

BeagleBuffer BeagleTreeLikelihood::borrowBuffer()
{
    assert(!availableBuffers_.empty());
    return BeagleBuffer(this);
}

void BeagleTreeLikelihood::returnBuffer(const int buffer, const bool check)
{
    assert(buffer < static_cast<int>(nBuffers_));
    auto it = bufferMap_.find(buffer);
    if(check && it == bufferMap_.end()) {
        throw std::runtime_error("Tried to return unknown buffer " + std::to_string(buffer));
    }
    if(it != bufferMap_.end()) {
        boost::clear_vertex(it->second, graph);
        boost::remove_vertex(it->second, graph);
        bufferMap_.erase(it);
        availableBuffers_.push(buffer);
    }
}

void BeagleTreeLikelihood::initialize(const bpp::SubstitutionModel& model,
                                      const bpp::DiscreteDistribution& rateDist,
                                      bpp::TreeTemplate<bpp::Node>& tree)
{
    this->rateDist_ = &rateDist;
    this->model_ = &model;
    this->tree_ = &tree;
    verifyInitialized();

    // Clear buffer maps
    distalNodeVertex_.clear();
    proxNodeVertex_.clear();
    bufferMap_.clear();
    leafVertex_.clear();
    graph = TGraph();

    std::vector<bool> isNonLeafBuffer(nBuffers_, true);
    for(const auto& p : leafBuffer_) {
        isNonLeafBuffer[p.second] = false;
    }

    assert(std::count(isNonLeafBuffer.begin(),
                      isNonLeafBuffer.end(),
                      false) ==
           static_cast<long>(leafBuffer_.size()));
    for(size_t i = 0; i < nBuffers_; i++) {
        if(isNonLeafBuffer[i]) {
            // The graph has been reset, so we don't need to remove the vertex.
            availableBuffers_.push(i);
        }
    }

    // Re-add leaf buffers
    for(const auto& p : leafBuffer_) {
        leafVertex_[p.first] = addBufferToGraph(VertexInfo{p.second, 0, false, true});
        bufferMap_[p.second] = leafVertex_[p.first];
    }
    assert(leafVertex_.size() == leafBuffer_.size() && "Leaf size does not match");

    loadRateDistribution(rateDist);
    loadSubstitutionModel(model);

    // Slide root to far right side of root branch
    tree.getRootNode()->getSon(0)->setDistanceToFather(tree.getRootNode()->getSon(0)->getDistanceToFather() +
                                                       tree.getRootNode()->getSon(1)->getDistanceToFather());
    tree.getRootNode()->getSon(1)->setDistanceToFather(0.0);

    // Fill buffer maps
    allocateDistalBuffers();
    allocateProximalBuffers();
    buildBufferDependencyGraph();
}

BeagleTreeLikelihood::TVertex BeagleTreeLikelihood::addBufferToGraph(const VertexInfo& info)
{
    if(bufferMap_.find(info.buffer) != bufferMap_.end())
        throw std::runtime_error("Buffer " + std::to_string(info.buffer) + " is already in graph.");
    TVertex vertex = boost::add_vertex(info, graph);
    bufferMap_[info.buffer] = vertex;
    return vertex;
}

void BeagleTreeLikelihood::addDependencies(const TVertex u,
                                           const TVertex v1, const double dist1,
                                           const TVertex v2, const double dist2)
{
    addDependency(u, v1, dist1);
    addDependency(u, v2, dist2);
}

bool BeagleTreeLikelihood::addDependency(const TVertex u, const TVertex v, const double dist)
{
    TEdge edge;
    bool added;
    if(u == v)
        throw std::runtime_error("Tried to add self dependency: " + std::to_string(u) + "->" + std::to_string(v));
    assert(u != v && "Cannot depend on itself!");

    std::tie(edge, added) = boost::add_edge(u, v, dist, graph);

    // Dependency already exists. overwrite current distance
    if(!added)
        graph[edge] = dist;

    return added;
}

void BeagleTreeLikelihood::allocateDistalBuffers()
{
    const std::vector<bpp::Node*> nodes = tree_->getNodes();
    for(const bpp::Node* n : nodes) {
        if(n->isLeaf()) {
            const std::string& name = n->getName();
            assert(leafVertex_.count(name) > 0);
            distalNodeVertex_[n] = leafVertex_.at(name);
        } else {
            assert(n->getNumberOfSons() == 2);
            assert(distalNodeVertex_.find(n) == distalNodeVertex_.end());
            distalNodeVertex_[n] = bufferMap_.at(getFreeBuffer());
        }
    }
}

void BeagleTreeLikelihood::allocateProximalBuffers()
{
    // Special handling at the root - proximal buffers for each child is the *distal* buffer for its sibling.
    const bpp::Node* root = tree_->getRootNode();
    for(size_t i = 0; i < root->getNumberOfSons(); i++) {
        const bpp::Node* n = root->getSon(i);
        proxNodeVertex_[n] = distalNodeVertex_.at(sibling(n));
    }
    for(const bpp::Node* n : tree_->getNodes()) {
        if(n->isLeaf() || n == root)
            continue;
        assert(n->getNumberOfSons() == 2);
        for(size_t i = 0; i < 2; ++i) {
            const bpp::Node* son = n->getSon(i);
            assert(distalNodeVertex_.find(son) != distalNodeVertex_.end());
            assert(proxNodeVertex_.find(son) == proxNodeVertex_.end());
            proxNodeVertex_[son] = bufferMap_.at(getFreeBuffer());
        }
    }
}

void BeagleTreeLikelihood::buildBufferDependencyGraph(bool allowExisting)
{
    if(!allowExisting && boost::num_edges(graph))
        throw new std::runtime_error("Expected an unconnected graph, got " +
                                     std::to_string(boost::num_edges(graph)) +
                                     " edges");
    // Distal
    // This is easy - each node just depends on the distal buffers of its children.
    for(const bpp::Node* n : tree_->getNodes()) {
        if(n->isLeaf())
            continue;
        assert(n->getNumberOfSons() == 2);
        for(size_t i = 0; i < n->getNumberOfSons(); i++) {
            const bpp::Node* son = n->getSon(i);

            const TVertex nodeV = distalNodeVertex_.at(n);
            const TVertex sonV = distalNodeVertex_.at(son);
            addDependency(nodeV, sonV, son->getDistanceToFather());
        }
    }

    // Proximal
    // Here, we traverse nodes, adding proximal nodes for the sons of each non-leaf, non-root node.
    for(const bpp::Node* parent : tree_->getNodes()) {
        if(parent->isLeaf() || parent == tree_->getRootNode())
            continue;
        assert(parent->getNumberOfSons() == 2);
        // The distal buffer for this node should already be calculated.
        assert(distalNodeVertex_.find(parent) != distalNodeVertex_.end());

        const TVertex parentVertex = proxNodeVertex_.at(parent);

        for(size_t i = 0; i < parent->getNumberOfSons(); ++i) {
            const bpp::Node* son = parent->getSon(i);
            const bpp::Node* sibling = siblings(son).at(0);

            assert(distalNodeVertex_.find(sibling) != distalNodeVertex_.end());
            assert(proxNodeVertex_.find(son) != proxNodeVertex_.end());
            const TVertex vertex = proxNodeVertex_.at(son);

            const TVertex siblingVertex = distalNodeVertex_.at(sibling);

            addDependency(vertex, siblingVertex, sibling->getDistanceToFather());
            double parentDist = parent->getDistanceToFather();
            if(parent->getFather() == tree_->getRootNode())
                parentDist += sts::online::sibling(parent)->getDistanceToFather();
            addDependency(vertex, parentVertex, parentDist);
        }
    }
}

size_t BeagleTreeLikelihood::registerLeaf(const bpp::Sequence& sequence)
{
    verifyInitialized();
    const int buffer = getFreeBuffer();
    graph[bufferMap_.at(buffer)].dirty = false;
    graph[bufferMap_.at(buffer)].leaf = true;
    if(leafVertex_.count(sequence.getName()) > 0)
        throw std::runtime_error("Duplicate sequence name: " + sequence.getName());

    leafBuffer_[sequence.getName()] = buffer;
    leafVertex_[sequence.getName()] = bufferMap_.at(buffer);
    const std::vector<double> seq_partials = get_partials(sequence, *model_, nRates_);
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

    nBeagleUpdateTransitionsCalls_ += operations.size();
    totalBeagleUpdateTransitionsCalls_ += operations.size();
}

void BeagleTreeLikelihood::updateTransitionsPartials(const TVertex vertex)
{
    // First, update dirty/clean status
    {
        std::unordered_map<TVertex, VertexInfo> scratch;
        BeagleMarkDirtyVisitor<TGraph> visitor(scratch);
        boost::vector_property_map<boost::default_color_type> colorVec(boost::num_vertices(graph));
        boost::depth_first_visit(graph, vertex, visitor, colorVec, visitor);
        visitor.update_graph(graph);
    }
    // Now update the partials of any dirty nodes
    {
        std::unordered_set<TVertex> visited;
        std::vector<BeagleOperation> operations;
        std::vector<int> nodeIndices;
        std::vector<double> branchLengths;
        BeagleUpdatePartialsVisitor<TGraph> visitor(operations, nodeIndices, branchLengths, visited);
        boost::vector_property_map<boost::default_color_type> colorVec(boost::num_vertices(graph));
        boost::depth_first_visit(graph, vertex, visitor, colorVec, visitor);
        visitor.update_graph(graph);
        assert(visited.size() == operations.size() && "Operation count did not match visited count");
        updateTransitionsPartials(operations, branchLengths, nodeIndices, BEAGLE_OP_NONE);
    }
}

std::vector<double> BeagleTreeLikelihood::calculateAttachmentLikelihood(const std::string& leafName,
                                                                        const bpp::Node* node,
                                                                        const double distalLength,
                                                                        const std::vector<double>& pendantBranchLengths)
{
    if(!leafVertex_.count(leafName))
        throw std::runtime_error("Unknown leaf: " + leafName);

    const TVertex leafVert = leafVertex_[leafName];
    const int leafBuf = graph[leafVert].buffer;

    const BeagleBuffer b = borrowBuffer();
    const TVertex vert = bufferMap_.at(b.value());

    double edgeLength = node->getDistanceToFather();
    if(node->getFather() == tree_->getRootNode())
        edgeLength += sibling(node)->getDistanceToFather();

    if(distalLength > edgeLength)
        throw std::runtime_error("Invalid distal length! " +
                                 std::to_string(distalLength) +
                                 " > " +
                                 std::to_string(edgeLength));

    const TVertex dist = distalNodeVertex_.at(node),
                  prox = proxNodeVertex_.at(node);

    addDependencies(vert, prox, distalLength, dist, edgeLength - distalLength);

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
    if(!leafVertex_.count(leafName))
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
    const TVertex vertex = distalNodeVertex_.at(node);
    updateTransitionsPartials(vertex);
    beagle_check(beagleGetPartials(beagleInstance_, graph[vertex].buffer, BEAGLE_OP_NONE, result.data()));
    return result;
}

LikelihoodVector BeagleTreeLikelihood::getProximalPartials(const bpp::Node* node)
{
    LikelihoodVector result(nRates_, nSites_, nStates_);
    const TVertex vertex = proxNodeVertex_.at(node);
    updateTransitionsPartials(vertex);
    beagle_check(beagleGetPartials(beagleInstance_, graph[vertex].buffer, BEAGLE_OP_NONE, result.data()));
    return result;
}

LikelihoodVector BeagleTreeLikelihood::getLeafPartials(const std::string& name) const
{
    const int buffer = graph[leafVertex_.at(name)].buffer;
    LikelihoodVector result(nRates_, nSites_, nStates_);
    beagle_check(beagleGetPartials(beagleInstance_, buffer, BEAGLE_OP_NONE, result.data()));
    return result;
}

int BeagleTreeLikelihood::getDistalBuffer(const bpp::Node* node) const
{
    return graph[distalNodeVertex_.at(node)].buffer;
}

int BeagleTreeLikelihood::getProximalBuffer(const bpp::Node* node) const
{
    return graph[proxNodeVertex_.at(node)].buffer;
}

int BeagleTreeLikelihood::getLeafBuffer(const std::string& name) const
{
    return graph[leafVertex_.at(name)].buffer;
}

void BeagleTreeLikelihood::invalidateAll()
{
    using TIterator = boost::graph_traits<TGraph>::vertex_iterator;
    TIterator it, end;
    for(std::tie(it, end) = boost::vertices(graph); it != end; ++it) {
        graph[*it].dirty = !graph[*it].leaf;
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
    TVertex root = distalNodeVertex_.at(tree_->getRootNode());
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
    for(const bpp::Node* n : postorder(tree_->getRootNode())) {
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

    const std::vector<bpp::Node*> nodes = tree_->getNodes();
    // Distal buffer
    for(const bpp::Node* n : nodes) {
        const int distalBuffer = distalNodeVertex_.at(n);
        out << "b" << distalBuffer << "[shape=none];\n";
        out << "b" << distalBuffer << " -> " << n->getId() << "[color=blue];\n";
    }
    const bpp::Node* root = tree_->getRootNode();
    for(size_t i = 0; i < root->getNumberOfSons(); i++) {
        const bpp::Node* n = root->getSon(i);
        const int proxBuffer = proxNodeVertex_.at(n);
        out << "b" << proxBuffer << "[shape=none];\n";
        out << "b" << proxBuffer << " -> " << n->getId() << "[color=red,style=dashed];\n";
    }
    for(const bpp::Node* n : nodes) {
        if(n->isLeaf() || n == root)
            continue;
        for(size_t i = 0; i < 2; ++i) {
            const bpp::Node* son = n->getSon(i);
            const int proxBuffer = proxNodeVertex_.at(son);
            out << "b" << proxBuffer << "[shape=none];\n";
            out << "b" << proxBuffer << " -> " << son->getId() << "[color=red,style=dashed];\n";
        }
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
    assert(buffer1 < static_cast<int>(nBuffers_) && buffer1 >= 0 && "Invalid buffer!");
    assert(buffer2 < static_cast<int>(nBuffers_) && buffer2 >= 0 && "Invalid buffer!");
    assert(freeBufferCount() >= 1);

    const BeagleBuffer b = borrowBuffer();
    const int scratchBuffer = b.value();
    assert(scratchBuffer != buffer1 && scratchBuffer != buffer2 &&
           "Reused buffer");
    addDependencies(bufferMap_.at(scratchBuffer),
                    bufferMap_.at(buffer1), 0,
                    bufferMap_.at(buffer2), d);

    updateTransitionsPartials(bufferMap_.at(scratchBuffer));

    return logLikelihood(scratchBuffer);
}

double BeagleTreeLikelihood::logLikelihood(const std::vector<double>& v)
{
    assert(v.size() == partialLength() && "unexpected partial length");
    const BeagleBuffer b = borrowBuffer();
    const int tmpBuffer = b.value();
    if(!(instanceDetails_.flags & BEAGLE_FLAG_SCALING_AUTO))
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
    const TVertex vertex = bufferMap_.at(buffer);
    std::vector<int> scaleBuffers;
    BeagleScaleFactorVisitor<TGraph> visitor(scaleBuffers);
    boost::vector_property_map<boost::default_color_type> colorVec(boost::num_vertices(graph));
    boost::depth_first_visit(graph, vertex, visitor, colorVec, visitor);

    beagle_check(beagleAccumulateScaleFactors(beagleInstance_,
                                              scaleBuffers.data(),
                                              scaleBuffers.size(),
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

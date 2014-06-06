/// \file beagle_tree_likelihood.h
/// \brief BEAGLE-Bio++ interface

#ifndef STS_LIKELIHOOD_BEAGLE_TREE_LIKELIHOOD_H
#define STS_LIKELIHOOD_BEAGLE_TREE_LIKELIHOOD_H

#include <Bpp/Phyl/TreeTemplate.h>

#include "libhmsbeagle/beagle.h"

#include <boost/graph/adjacency_list.hpp>

#include <iostream>
#include <stack>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Forwards
namespace bpp {
class SiteContainer;
class SubstitutionModel;
class Sequence;
class DiscreteDistribution;
}

namespace sts { namespace online {

class LikelihoodVector;
class BeagleBuffer;

/// \brief Beagle-Bio++ interface
///
/// Calculates the likelihood of a Bio++ tree using BEAGLE, with dependency tracking between nodes,
/// so after a branch length change, only a subset of nodes are re-peeled.
///
/// \todo Invalidate on model / rate distribution changes
///
/// \todo Nodes are hashed using a combination of branch length and child addresses. Is that sufficient?
///
/// Currently, constructed with an alignment *containing all sequences that will be used, including sequences / to be
/// added*, a substitution model, and a rate distribution, see #BeagleTreeLikelihood.
///
/// Prior to use, #initialize must be called with the current tree, substitution model, and rate distribution.
///
/// \note During initialization, the input tree is changed such that the root is as far to the right as possible on the
/// current branch, so for the topology:
/// <pre>
///       root
///  l1   /  \  l2
///      n1   n1
/// </pre>
/// Initialization will set <c>l1 = l1 + l2, l2 = 0</c>.
class BeagleTreeLikelihood
{
    friend class BeagleBuffer;
public:
    /// \brief Constructor
    ///
    /// Initializes a BEAGLE instance, fills partials vectors for every sequence in \c sites
    /// \param sites Alignment. May contain more sequences than taxa present in future calls to #initialize.
    /// Sufficient buffers are allocated for a fully-specified tree.
    /// \param model Default subsitution model. Used solely for filling partials vector. A model with the same number of
    /// states should later be specified via #loadSubstitutionModel.
    /// \param rateDist Discrete rate distribution. Used solely for filling partials vector. The actual rate
    /// distribution associated with a tree should be specified via #loadRateDistribution.
    /// \param extraBufferCount Number of spare buffers to allocate.
    BeagleTreeLikelihood(const bpp::SiteContainer& sites,
                         const bpp::SubstitutionModel& model,
                         const bpp::DiscreteDistribution& rateDist,
                         const size_t extraBufferCount=3);

    /// \brief Move constructor
    BeagleTreeLikelihood(BeagleTreeLikelihood&& other) = default;
    BeagleTreeLikelihood(const BeagleTreeLikelihood& other) = delete;

    /// Destructor - frees BEAGLE resources
    virtual ~BeagleTreeLikelihood();

    /// \brief Initialize the beagle_instance for a model, rate distribution, and tree
    void initialize(const bpp::SubstitutionModel& model,
                    const bpp::DiscreteDistribution& rateDist,
                    bpp::TreeTemplate<bpp::Node>& tree);

    /// \brief Calculate the likelihood of a tree.
    ///
    /// \returns log-likelihood.
    double calculateLogLikelihood();

    /// \brief Gets the BEAGLE instance ID associated with this instance.
    int beagleInstance() const { return beagleInstance_; };
    /// \brief Number of buffers allocated
    size_t numberOfBuffers() const { return nBuffers_; };
    /// \brief Length of a single partial likelihood vector
    size_t partialLength() const { return nSites_ * nStates_ * nRates_; };

    /// Get the partials for the distal side of an edge
    LikelihoodVector getDistalPartials(const bpp::Node* node);
    /// Get the partials for the proximal side of an edge
    LikelihoodVector getProximalPartials(const bpp::Node* node);

    /// Get the BEAGLE buffer index for the distal side of an edge
    int getDistalBuffer(const bpp::Node* node) const;
    /// Get the BEAGLE buffer index for the proximal side of an edge
    int getProximalBuffer(const bpp::Node* node) const;
    /// Get the mid-edge buffer associated with a Node
    int getMidEdgeBuffer(const bpp::Node* node) const;
    /// Get the BEAGLE buffer index a leaf by name
    int getLeafBuffer(const std::string& name) const;

    /// Borrow a buffer from the set of free buffers.
    BeagleBuffer borrowBuffer();

    size_t freeBufferCount() const { return availableBuffers.size(); };

    typedef std::pair<const bpp::Node*, int> NodePartials;

    /// \brief Get a partial vector for the middle of each edge
    std::vector<NodePartials> getMidEdgePartials();

    /// \brief Get the partials for a given taxon
    ///
    /// \param name Taxon name
    LikelihoodVector getLeafPartials(const std::string& name) const;

    /// \brief pair representing a node and a distance from the proximal side of the edge for attachment
    typedef std::pair<bpp::Node*, double> AttachmentLocation;

    /// \brief Calculate the attachment likelihood of the leaf #leafName at various locations.
    ///
    /// \param leafName Name of the leaf - must be registered.
    /// \param node Node representing the edge on which to attach leafName
    /// \param distalLength distance from pendant side of edge
    /// \param pendant Pendant edge length
    ///
    /// \return Log-likelihood associated with the attachment
    std::vector<double> calculateAttachmentLikelihood(const std::string& leafName,
                                                      const bpp::Node* node,
                                                      const double distalLength,
                                                      const std::vector<double>& pendantBranchLengths = std::vector<double>{0.0});

    /// \brief Calculate the attachment likelihood of the leaf #leafName at various locations.
    ///
    /// \param leafName Name of the leaf - must be registered.
    /// \param attachmentLocations Pairs consisting of a node and distance from proximal side.
    /// \return Log likelihood associated with each attachment location.
    std::vector<std::vector<double>> calculateAttachmentLikelihoods(const std::string& leafName,
                                                                    const std::vector<AttachmentLocation>& attachmentLocations,
                                                                    const std::vector<double> pendantBranchLengths = std::vector<double>{0.0});

    /// \brief Invalidate a single node (indicating that that node, and parents should be re-peeled).
    void invalidate(const bpp::Node* node);
    /// \brief Invalidate all nodes. Full re-peel will be performed on next likelihood call.
    void invalidateAll();

    /// \brief Log-likelihood of two partials vectors connected by a branch of length \c d
    double logDot(const std::vector<double>& v1, const std::vector<double>& v2, const double d = 0.0);
    /// \brief Log-likelihood of two partials vectors connected by a branch of length \c d
    double logDot(const std::vector<double>& v, const int buffer, const double d = 0.0);
    /// \brief Log-likelihood of two partials vectors connected by a branch of length \c d
    double logDot(const int buffer1, const int buffer2, const double d = 0.0);

    /// Calculate the summed log-likelihood of a partials vector
    double logLikelihood(const std::vector<double>& v);
    /// Calculate the summed log-likelihood of a partials buffer
    double logLikelihood(const int buffer);


protected:
    /// \brief Load eigendecomposition of \c model
    ///
    /// \param model Model, with same number of states as specified in the constructor
    void loadSubstitutionModel(const bpp::SubstitutionModel& model);

    /// \brief Load rates and weights from \c rate_dist

    /// \param rate_dist Rate distribution, with same number of categories as passed in constructor.
    void loadRateDistribution(const bpp::DiscreteDistribution& rate_dist);

    /// \brief Calculate distal partial vectors for every internal node in the tree.
    void calculateDistalPartials();

    /// \brief Calculate proximal partial vectors for every internal node in the tree.
    void calculateProximalPartials();

    /// \brief Update the transition matrices and partial likelihood vectors.
    ///
    /// \param operations A vector of beagle operations
    /// \param branchLengths a vector of length <c>2 * operations.size()</c>, containing branch lengths for the two
    /// children of each node in \c operations.
    /// \param nodeIndices Buffer indices for the nodes referred to in \c branch_lengths.
    /// \param rootBuffer Root buffer - used for scaling
    void updateTransitionsPartials(const std::vector<BeagleOperation>& operations,
                                   const std::vector<double>& branchLengths,
                                   const std::vector<int>& nodeIndices,
                                   const int rootBuffer);

    void accumulateScaleFactors(const std::vector<BeagleOperation>& operations, const int scalingBuffer);

    /// \brief Register a leaf sequence
    size_t registerLeaf(const bpp::Sequence& sequence);

    /// Borrow a buffer from the set of free buffers.
    /// Will not be re-used until #returnBuffer is called.
    int getFreeBuffer();
    /// Return a buffer which is no longer needed.
    ///
    /// \param buffer Buffer index
    /// \param check should function assert that \c buffer has been allocated via #getFreeBuffer?
    void returnBuffer(const int buffer, const bool check=true);

private:
    // Buffer dependencies
    struct VertexInfo
    {
        int buffer;
        size_t hash;
        bool dirty;
    };
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, VertexInfo, double> TGraph;
    typedef boost::graph_traits<TGraph>::vertex_descriptor TVertex;
    /// Dependency graph between buffers
    TGraph graph;

    void verifyInitialized() const;
    int beagleInstance_;

    const size_t nSites_;
    const size_t nStates_;
    const size_t nRates_;
    const size_t nSeqs_;
    const size_t nBuffers_;

    // Buffer tracking
    std::stack<int> availableBuffers;
    std::unordered_set<int> usedBuffers;

    BeagleInstanceDetails instanceDetails;

    /// Map from leaf name to BEAGLE buffer
    std::unordered_map<std::string, TVertex> leafBuffer;

    /// Model stuff
    bpp::DiscreteDistribution const* rateDist;
    bpp::SubstitutionModel const* model;
    bpp::TreeTemplate<bpp::Node>* tree;

    /// Map from node to the BEAGLE buffer for its distal partial vector
    std::unordered_map<const bpp::Node*, TVertex> distalNodeBuffer;
    /// Map from node to the BEAGLE buffer for its proximal partial vector
    std::unordered_map<const bpp::Node*, TVertex> proxNodeBuffer;
    /// Map from node to the BEAGLE buffer for the middle of the edge above the node
    std::unordered_map<const bpp::Node*, TVertex> midEdgeNodeBuffer;

    std::unordered_map<int, TVertex> bufferMap;

    /// Map from a node to a hash of its state last time its distal likelihood vector was calculated
    std::unordered_map<const bpp::Node*, size_t> distalNodeState;
    /// Map from a node to a hash of its state last time its proximal likelihood vector was calculated
    std::unordered_map<const bpp::Node*, size_t> proxNodeState;

    /// For testing, mostly. Writes a graph with node numbers, prox / distal buffer indices.
    void toDot(std::ostream& out) const;

    TVertex addBufferToGraph(const VertexInfo& info);
    void allocateDistalBuffers();
    void allocateProximalBuffers();
    void allocateMidEdgeBuffers();
    void buildBufferDependencyGraph();
};

/// Representation of a Beagle Buffer.
/// Buffer is returned when destroyed.
class BeagleBuffer
{
public:
    BeagleBuffer(BeagleTreeLikelihood* btl) :
        instance_(btl),
        value_(instance_->getFreeBuffer()) {};
    BeagleBuffer(const BeagleBuffer&) = delete;
    BeagleBuffer& operator=(const BeagleBuffer&) = delete;
    BeagleBuffer(BeagleBuffer&& other) = default;
    BeagleBuffer& operator=(BeagleBuffer&& other) = default;
    ~BeagleBuffer() { instance_->returnBuffer(value_); };

    /// Gets the buffer index
    int value() const { return value_; };
private:
    BeagleTreeLikelihood* instance_;
    int value_;
};

}}

#endif

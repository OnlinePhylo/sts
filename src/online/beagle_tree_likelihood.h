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
#include <vector>

#include "attachment_location.h"

#include <Bpp/Phyl/SitePatterns.h>

// Forwards
namespace bpp {
class SiteContainer;
class SitePatterns;
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
/// \todo Notify on branch length change?
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
///
/// All computations should be lazy - BeagleTreeLikelihood instances store a graph managing dependencies between buffers
/// in BEAGLE.
/// \note Node dependencies are calculated in BeagleTreeLikelihood::initialize.
/// Changing the tree *after* initializing is not supported (though it shouldn't be too hard - see
/// BeagleTreeLikelihood::buildBufferDependencyGraph).
class BeagleTreeLikelihood
{
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
    BeagleTreeLikelihood(const bpp::SitePatterns& patterns,
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
    size_t partialLength() const { return patternCount_ * nStates_ * nRates_; };

    /// \brief Total number of operations passed to beagleUpdateTransitionMatrices
    size_t numberOfBeagleUpdateTransitionsCalls() const { return nBeagleUpdateTransitionsCalls_; }

    static size_t totalBeagleUpdateTransitionsCalls() { return totalBeagleUpdateTransitionsCalls_; }

    /// \brief Get the partials for the distal side of an edge
    LikelihoodVector getDistalPartials(const bpp::Node* node);
    /// \brief Get the partials for the proximal side of an edge
    LikelihoodVector getProximalPartials(const bpp::Node* node);

    /// \brief Get the BEAGLE buffer index for the distal side of an edge
    int getDistalBuffer(const bpp::Node* node) const;
    /// \brief Get the BEAGLE buffer index for the proximal side of an edge
    int getProximalBuffer(const bpp::Node* node) const;
    /// \brief Get the BEAGLE buffer index a leaf by name
    int getLeafBuffer(const std::string& name) const;

    /// \brief Borrow a buffer from the set of free buffers.
    /// \pre BeagleTreeLikelihood::freeBufferCount() > 0
    std::unique_ptr<BeagleBuffer> borrowBuffer();

    size_t freeBufferCount() const { return availableBuffers_.size(); };

    using NodePartials = std::pair<const bpp::Node*, int>;

    /// \brief Get the partials for a given taxon
    ///
    /// \param name Taxon name
    LikelihoodVector getLeafPartials(const std::string& name) const;

    /// \brief Calculate the attachment likelihood of the leaf #leafName at various locations.
    ///
    /// \param leafName Name of the leaf - must be registered.
    /// \param node Node representing the edge on which to attach leafName
    /// \param distalLength distance from *root* side of edge
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
    friend class BeagleBuffer;

    // typedefs for dependency tracking
    // Description of a vertex.
    // Each vertex is associated with a BEAGLE buffer.
    struct VertexInfo
    {
        int buffer;
        size_t hash;
        bool dirty;
        bool leaf;
    };

    /// Dependency graph between buffers. Edges are described by a double indicating length,
    /// vertices are described by a VertexInfo
    /// Vertices are stored as lists - this prevents invalidating TVertex's on
    /// vertex removal.
    /// Edges are stored as vectors, so edge iterators *are* invalidated on edge removal.
    using TGraph = boost::adjacency_list<boost::listS, boost::vecS, boost::directedS, VertexInfo, double>;
    using TVertex = boost::graph_traits<TGraph>::vertex_descriptor;
    using TEdge = boost::graph_traits<TGraph>::edge_descriptor;

    /// For testing, mostly. Writes a graph with node numbers, prox / distal buffer indices.
    void toDot(std::ostream& out) const;

    /// \brief Add a dependency of \c u on \c v1 and \c v2, with associated distances \c dist1 and \c dist2.
    ///
    /// If u is already dependent on \c v1 or \c v2, the edge is updated to have the new distance.
    TVertex addBufferToGraph(const VertexInfo& info);

    /// \brief introduce dependencies of `u` on `v1` and `v2`
    void addDependencies(const TVertex u,
                         const TVertex v1, const double dist1,
                         const TVertex v2, const double dist2);

    /// \brief Add a dependency of `u` on `v` with given distance.
    ///
    /// If `u` already depends on `v`, the distance is updated and no other changes are made.
    /// \return whether or not a new edge was introduced (`false` indicates edge update)
    bool addDependency(const TVertex u, const TVertex v, const double dist);
    void allocateDistalBuffers();
    void allocateProximalBuffers();

    /// \brief Add dependencies between buffers based on the input tree.
    void buildBufferDependencyGraph(bool allowExisting = false);

    /// \brief Update transitions matrices and partials vectors of *all dirty nodes* in the tree rooted at `vertex`
    void updateTransitionsPartials(const TVertex vertex);

    void verifyInitialized() const;

    // Buffer dependency tracking
    /// Dependency graph between buffers
    TGraph graph;

    int beagleInstance_;

    const size_t patternCount_;
    const size_t nStates_;
    const size_t nRates_;
    size_t nSeqs_;
    size_t nBuffers_;
    size_t nBeagleUpdateTransitionsCalls_;
    static size_t totalBeagleUpdateTransitionsCalls_;

    // Buffer tracking
    std::stack<int> availableBuffers_;

    BeagleInstanceDetails instanceDetails_;

    /// Map from leaf name to leaf vertex in graph
    std::unordered_map<std::string, TVertex> leafVertex_;
    /// Map from leaf name to leaf buffer in BEAGLE
    std::unordered_map<std::string, int> leafBuffer_;

    /// Model stuff
    const bpp::SitePatterns& _patterns;
    bpp::DiscreteDistribution const* rateDist_;
    bpp::SubstitutionModel const* model_;
    bpp::TreeTemplate<bpp::Node>* tree_;

    /// Map from node to the BEAGLE buffer for its distal partial vector
    std::unordered_map<const bpp::Node*, TVertex> distalNodeVertex_;
    /// Map from node to the BEAGLE buffer for its proximal partial vector
    std::unordered_map<const bpp::Node*, TVertex> proxNodeVertex_;

    /// Map from a buffer to a vertex in `graph`
    std::unordered_map<int, TVertex> bufferMap_;
};

/// Representation of a Beagle Buffer.
/// Buffer is returned when destroyed.
class BeagleBuffer
{
public:
    BeagleBuffer(BeagleTreeLikelihood* btl) :
        instance_(btl),
        value_(instance_->getFreeBuffer()) { };
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

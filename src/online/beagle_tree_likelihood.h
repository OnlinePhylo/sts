/// \file beagle_tree_likelihood.h
/// \brief BEAGLE-Bio++ interface

#ifndef STS_LIKELIHOOD_BEAGLE_TREE_LIKELIHOOD_H
#define STS_LIKELIHOOD_BEAGLE_TREE_LIKELIHOOD_H

#include <Bpp/Phyl/TreeTemplate.h>

#include "libhmsbeagle/beagle.h"

#include <iostream>
#include <stack>
#include <string>
#include <unordered_map>
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
                         const size_t extraBufferCount=2);

    /// \brief Move constructor
    BeagleTreeLikelihood(BeagleTreeLikelihood&& other);

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
    int getDistalBuffer(const bpp::Node* node);
    /// Get the BEAGLE buffer index for the proximal side of an edge
    int getProximalBuffer(const bpp::Node* node);
    /// Get the BEAGLE buffer index a leaf by name
    int getLeafBuffer(const std::string& name);

    typedef std::pair<const bpp::Node*, LikelihoodVector> NodePartials;

    /// \brief Get a partial vector for the middle of each edge
    std::vector<NodePartials> getMidEdgePartials();

    /// \brief Get the partials for a given taxon
    ///
    /// \param name Taxon name
    LikelihoodVector getLeafPartials(const std::string& name);

    /// \brief Invalidate a single node (indicating that that node, and parents should be re-peeled).
    void invalidate(const bpp::Node* node);
    /// \brief Invalidate all nodes. Full re-peel will be performed on next likelihood call.
    void invalidateAll();

    inline const std::vector<int>& getScratchBuffers() const { return scratchBuffers_; }
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
private:
    void verifyInitialized() const;
    int beagleInstance_;

    const size_t nSites_;
    const size_t nStates_;
    const size_t nRates_;
    const size_t nSeqs_;
    const size_t nBuffers_;

    const size_t nScratchBuffers_;
    std::vector<int> scratchBuffers_;

    BeagleInstanceDetails instanceDetails;

    /// Map from leaf name to BEAGLE buffer
    std::unordered_map<std::string, int> leafBuffer;

    /// Model stuff
    bpp::DiscreteDistribution const* rateDist;
    bpp::SubstitutionModel const* model;
    bpp::TreeTemplate<bpp::Node>* tree;

    /// Map from node to the BEAGLE buffer for its distal partial vector
    std::unordered_map<const bpp::Node*, int> distalNodeBuffer;
    /// Map from node to the BEAGLE buffer for its proximal partial vector
    std::unordered_map<const bpp::Node*, int> proxNodeBuffer;

    /// Map from a node to a hash of its state last time its distal likelihood vector was calculated
    std::unordered_map<const bpp::Node*, size_t> distalNodeState;
    /// Map from a node to a hash of its state last time its proximal likelihood vector was calculated
    std::unordered_map<const bpp::Node*, size_t> proxNodeState;

    /// For testing, mostly. Writes a graph with node numbers, prox / distal buffer indices.
    void toDot(std::ostream& out) const;
};

}}

#endif

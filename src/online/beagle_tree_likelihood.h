/// \file beagle_tree_likelihood.h
/// \brief BEAGLE-Bio++ interface

#ifndef STS_LIKELIHOOD_BEAGLE_TREE_LIKELIHOOD_H
#define STS_LIKELIHOOD_BEAGLE_TREE_LIKELIHOOD_H

#include <Bpp/Phyl/TreeTemplate.h>

#include "libhmsbeagle/beagle.h"

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

class Likelihood_vector;

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
/// added*, a substitution model, and a rate distribution, see #Beagle_tree_likelihood.
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
class Beagle_tree_likelihood
{

public:
    /// \brief Constructor
    ///
    /// Initializes a BEAGLE instance, fills partials vectors for every sequence in \c sites
    /// \param sites Alignment. May contain more sequences than taxa present in future calls to #load_tree.
    /// Sufficient buffers are allocated for a fully-specified tree.
    /// \param model Default subsitution model. Used solely for filling partials vector. A model with the same number of
    /// states should later be specified via #load_substitution_model.
    /// \param rate_distribution Discrete rate distribution. Used solely for filling partials vector. The actual rate
    /// distribution associated with a tree should be specified via #load_rate_distribution.
    Beagle_tree_likelihood(const bpp::SiteContainer& sites,
                           const bpp::SubstitutionModel& model,
                           const bpp::DiscreteDistribution& rate_dist,
                           const size_t extra_buffer_count=2);

    /// \brief Move constructor
    Beagle_tree_likelihood(Beagle_tree_likelihood&& other);

    /// Destructor - frees BEAGLE resources
    virtual ~Beagle_tree_likelihood();

    /// \brief Initialize the beagle_instance for a model, rate distribution, and tree
    void initialize(const bpp::SubstitutionModel& model,
                    const bpp::DiscreteDistribution& rate_dist,
                    bpp::TreeTemplate<bpp::Node>& tree);

    /// \brief Calculate the likelihood of a tree.
    ///
    /// \param tree Tree, containing taxa matching sequences given in the constructor.
    /// \returns log-likelihood.
    double calculate_log_likelihood();

    /// \brief Gets the BEAGLE instance ID associated with this instance.
    int get_beagle_instance() const { return beagle_instance; };
    /// \brief Number of buffers allocated
    size_t get_n_buffers() const { return n_buffers; };
    /// \brief Length of a single partial likelihood vector
    size_t get_partial_length() const { return n_sites * n_states * n_rates; };

    /// Get the partials for the distal side of an edge
    Likelihood_vector get_distal_partials(const bpp::Node* node);
    /// Get the partials for the proximal side of an edge
    Likelihood_vector get_proximal_partials(const bpp::Node* node);

    /// Get the BEAGLE buffer index for the distal side of an edge
    int get_distal_buffer(const bpp::Node* node);
    /// Get the BEAGLE buffer index for the proximal side of an edge
    int get_proximal_buffer(const bpp::Node* node);
    /// Get the BEAGLE buffer index a leaf by name
    int get_leaf_buffer(const std::string& name);

    typedef std::pair<const bpp::Node*, Likelihood_vector> Node_partials;

    /// \brief Get a partial vector for the middle of each edge
    std::vector<Node_partials> get_mid_edge_partials();

    /// \brief Get the partials for a given taxon
    ///
    /// \param name Taxon name
    Likelihood_vector get_leaf_partials(const std::string& name);

    /// \brief Invalidate a single node (indicating that that node, and parents should be re-peeled).
    void invalidate(const bpp::Node* node);
    /// \brief Invalidate all nodes. Full re-peel will be performed on next likelihood call.
    void invalidate_all();

    inline const std::vector<int>& get_scratch_buffers() const { return scratch_buffers; }
protected:
    /// \brief Load eigendecomposition of \c model
    ///
    /// \param model Model, with same number of states as specified in the constructor
    void load_substitution_model(const bpp::SubstitutionModel& model);

    /// \brief Load rates and weights from \c rate_dist

    /// \param rate_dist Rate distribution, with same number of categories as passed in constructor.
    void load_rate_distribution(const bpp::DiscreteDistribution& rate_dist);

    /// \brief Calculate distal partial vectors for every internal node in the tree.
    void calculate_distal_partials();

    /// \brief Calculate proximal partial vectors for every internal node in the tree.
    void calculate_proximal_partials();

    /// \brief Update the transition matrices and partial likelihood vectors.
    ///
    /// \param operations A vector of beagle operations
    /// \param branch_lengths a vector of length <c>2 * operations.size()</c>, containing branch lengths for the two
    /// children of each node in \c operations.
    /// \param node_indices Buffer indices for the nodes referred to in \c branch_lengths.
    void update_transitions_partials(const std::vector<BeagleOperation>& operations,
                                     const std::vector<double>& branch_lengths,
                                     const std::vector<int>& node_indices,
                                     const int root_buffer);
    void accumulate_scale_factors(const std::vector<BeagleOperation>& operations, const int scaling_buffer);

    /// \brief Register a leaf sequence
    size_t register_leaf(const bpp::Sequence& sequence);
private:
    void verify_initialized() const;
    int beagle_instance;

    const size_t n_sites;
    const size_t n_states;
    const size_t n_rates;
    const size_t n_seqs;
    const size_t n_buffers;

    const size_t n_scratch_buffers;
    std::vector<int> scratch_buffers;

    BeagleInstanceDetails instance_details;

    /// Map from leaf name to BEAGLE buffer
    std::unordered_map<std::string, int> leaf_buffer;

    /// Model stuff
    bpp::DiscreteDistribution const* rate_dist;
    bpp::SubstitutionModel const* model;
    bpp::TreeTemplate<bpp::Node>* tree;

    /// Map from node to the BEAGLE buffer for its distal partial vector
    std::unordered_map<const bpp::Node*, int> distal_node_buffer;
    /// Map from node to the BEAGLE buffer for its proximal partial vector
    std::unordered_map<const bpp::Node*, int> prox_node_buffer;

    /// Map from a node to a hash of its state last time its distal likelihood vector was calculated
    std::unordered_map<const bpp::Node*, size_t> distal_node_state;
    /// Map from a node to a hash of its state last time its proximal likelihood vector was calculated
    std::unordered_map<const bpp::Node*, size_t> prox_node_state;
};

}}

#endif

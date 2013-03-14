/// \file beagle_tree_likelihood.h
/// \brief BEAGLE-Bio++ interface

#ifndef STS_LIKELIHOOD_BEAGLE_TREE_LIKELIHOOD_H
#define STS_LIKELIHOOD_BEAGLE_TREE_LIKELIHOOD_H

#include <Bpp/Phyl/TreeTemplate.h>

#include "libhmsbeagle/beagle.h"

#include <stack>
#include <string>
#include <unordered_map>

// Forwards
namespace bpp {
class SiteContainer;
class SubstitutionModel;
class Sequence;
class DiscreteDistribution;
}

namespace sts { namespace online {

/// \brief Beagle-Bio++ interface
///
/// Calculates the likelihood of a Bio++ tree.
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
                           const bpp::DiscreteDistribution& rate_dist);

    /// Destructor - frees BEAGLE resources
    virtual ~Beagle_tree_likelihood();

    /// \brief Load eigendecomposition of \c model
    ///
    /// \param model Model, with same number of states as specified in the constructor
    void load_substitution_model(const bpp::SubstitutionModel& model);

    /// \brief Load rates and weights from \c rate_dist

    /// \param rate_dist Rate distribution, with same number of categories as passed in constructor.
    void load_rate_distribution(const bpp::DiscreteDistribution& rate_dist);

    /// \brief Calculate the likelihood of a tree.
    ///
    /// \param tree Tree, containing taxa matching sequences given in the constructor.
    /// \returns log-likelihood.
    double calculate_log_likelihood(const bpp::TreeTemplate<bpp::Node>& tree);

    /// \brief Gets the BEAGLE instance ID associated with this instance.
    int get_beagle_instance() const { return beagle_instance; };
private:
    size_t register_leaf(const bpp::Sequence& sequence, const bpp::SubstitutionModel& model);
    void verify_initialized() const;
    int beagle_instance;
    const size_t n_sites;
    const size_t n_states;
    const size_t n_rates;
    const size_t n_seqs;
    const size_t n_buffers;
    BeagleInstanceDetails instance_details;
    std::unordered_map<std::string, int> leaf_buffer;
};

}}

#endif

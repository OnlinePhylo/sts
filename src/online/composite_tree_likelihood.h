#ifndef STS_ONLINE_COMPOSITE_TREE_LIKELIHOOD_H
#define STS_ONLINE_COMPOSITE_TREE_LIKELIHOOD_H

#include <Bpp/Phyl/TreeTemplate.h>

#include <functional>
#include <memory>
#include <vector>

namespace bpp {
class DiscreteDistribution;
class Node;
class SiteContainer;
class SubstitutionModel;
}

namespace sts { namespace online {

// Forwards
class AttachmentLikelihood;
class BeagleTreeLikelihood;
class TripodOptimizer;

/// Log-likelihood function for a tree
typedef std::function<double(bpp::TreeTemplate<bpp::Node>&)> TreeLogLikelihood;

/// Composite tree likelihood, consisting of a #sts::online::BeagleTreeLikelihood, and additional likelihoods (e.g.,
/// priors).
/// #sts::online::BeagleTreeLikelihood is treated differently due to additional functionality available.
class CompositeTreeLikelihood
{
public:
    explicit CompositeTreeLikelihood(std::shared_ptr<BeagleTreeLikelihood> calculator);
    CompositeTreeLikelihood(std::shared_ptr<BeagleTreeLikelihood> calculator,
                            std::vector<TreeLogLikelihood> additionalLogLikes);

    /// Add a tree likelihood function
    void add(TreeLogLikelihood like);

    /// \brief Initialize for a {model, rate_dist, tree}
    ///
    /// <b>Must be called before #CompositeTreeLikelihood::operator()()</b>
    void initialize(const bpp::SubstitutionModel& model,
                    const bpp::DiscreteDistribution& rate_dist,
                    bpp::TreeTemplate<bpp::Node>& tree);

    /// Calculate the sum of log-likelihoods. Alias for #logLikelihood
    double operator()();

    /// Calculate the sum of log likelihoods
    double logLikelihood();

    const std::vector<double> edgeLogLikelihoods(const std::string& leaf_name, const std::vector<double>& pendant_lengths);
    TripodOptimizer createOptimizer(const bpp::Node* insertEdge, const std::string& newLeafName);

private:
    friend class sts::online::AttachmentLikelihood;

    std::shared_ptr<BeagleTreeLikelihood> calculator_;
    std::vector<TreeLogLikelihood> additionalLogLikes;

    bpp::TreeTemplate<bpp::Node>* tree;
};

}} // Namespace

#endif

#ifndef STS_ONLINE_COMPOSITE_TREE_LIKELIHOOD_H
#define STS_ONLINE_COMPOSITE_TREE_LIKELIHOOD_H

#include <Bpp/Phyl/TreeTemplate.h>

#include <functional>
#include <memory>
#include <vector>

namespace bpp {
class SiteContainer;
class SubstitutionModel;
class DiscreteDistribution;
}

namespace sts { namespace online {

// Forwards
class BeagleTreeLikelihood;

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

    /// Gets the BEAGLE likelihood calculator
    inline std::shared_ptr<BeagleTreeLikelihood> calculator() { return calculator_; }

private:
    std::shared_ptr<BeagleTreeLikelihood> calculator_;
    std::vector<TreeLogLikelihood> additionalLogLikes;

    bpp::TreeTemplate<bpp::Node>* tree;
};

}} // Namespace

#endif

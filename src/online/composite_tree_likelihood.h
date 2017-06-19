#ifndef STS_ONLINE_COMPOSITE_TREE_LIKELIHOOD_H
#define STS_ONLINE_COMPOSITE_TREE_LIKELIHOOD_H

#include <Bpp/Phyl/TreeTemplate.h>

#include <functional>
#include <memory>
#include <vector>

#include "flexible_tree_likelihood.h"
#include "prior.h"

namespace bpp {
class DiscreteDistribution;
class Node;
class SiteContainer;
class SubstitutionModel;
}

namespace sts { namespace online {

// Forwards
class AttachmentLikelihood;
class TripodOptimizer;

/// Log-likelihood function for a tree
typedef std::function<double(bpp::TreeTemplate<bpp::Node>&)> TreeLogLikelihood;

/// Composite tree likelihood, consisting of a #sts::online::BeagleTreeLikelihood, and additional likelihoods (e.g.,
/// priors).
/// #sts::online::BeagleTreeLikelihood is treated differently due to additional functionality available.
class CompositeTreeLikelihood
{
public:
    explicit CompositeTreeLikelihood(std::shared_ptr<FlexibleTreeLikelihood> calculator);
    CompositeTreeLikelihood(std::shared_ptr<FlexibleTreeLikelihood> calculator,
                            std::vector<TreeLogLikelihood> additionalLogLikes);

    /// Add a tree likelihood function
    void add(TreeLogLikelihood like);
    
    void add(std::unique_ptr<Prior>& prior);

    /// Sum additional log-likelihood functions
    double sumAdditionalLogLikes() const;

    /// \brief Initialize for a {model, rate_dist, tree}
    ///
    /// <b>Must be called before #CompositeTreeLikelihood::operator()()</b>
    void initialize(const bpp::SubstitutionModel& model,
                    const bpp::DiscreteDistribution& rate_dist,
                    bpp::TreeTemplate<bpp::Node>& tree);

    /// Calculate the sum of log-likelihoods. Alias for #logLikelihood
    double operator()();
    
    double operator()(const bpp::Node& distal, std::string taxonName, double pendantLength, double distalLength, double proximalLength);

    /// Calculate the sum of log likelihoods
    double logLikelihood();
    
    double logLikelihood(const bpp::Node& distal, std::string taxonName, double pendantLength, double distalLength, double proximalLength);
    
    void calculatePendantDerivatives(const bpp::Node& distal, std::string taxonName, double pendantLength, double distalLength, double proximalLength, double* d1, double* d2);
    
    void calculateDistalDerivatives(const bpp::Node& distal, std::string taxonName, double pendantLength, double distalLength, double proximalLength, double* d1, double* d2);

private:
    friend class sts::online::AttachmentLikelihood;

    std::shared_ptr<FlexibleTreeLikelihood> calculator_;
    std::vector<TreeLogLikelihood> additionalLogLikes_;
    
    std::vector<std::unique_ptr<Prior>> _priors;

    bpp::TreeTemplate<bpp::Node>* tree_;
};

}} // Namespace

#endif

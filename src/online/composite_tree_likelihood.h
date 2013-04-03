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
class Beagle_tree_likelihood;

/// Log-likelihood function for a tree
typedef std::function<double(bpp::TreeTemplate<bpp::Node>&)> Tree_log_likelihood;

/// Composite tree likelihood, consisting of a #sts::online::Beagle_tree_likelihood, and additional likelihoods (e.g., priors).
/// #sts::online::Beagle_tree_likelihood is treated differently due to additional functionality available.
class Composite_tree_likelihood
{
public:
    explicit Composite_tree_likelihood(std::shared_ptr<Beagle_tree_likelihood> calculator);
    Composite_tree_likelihood(std::shared_ptr<Beagle_tree_likelihood> calculator,
                              std::vector<Tree_log_likelihood> additional_log_likes);

    /// Add a tree likelihood function
    void add(Tree_log_likelihood like);

    /// \brief Initialize for a {model, rate_dist, tree}
    ///
    /// <b>Must be called before #Composite_tree_likelihood::operator()</b>
    void initialize(const bpp::SubstitutionModel& model,
                    const bpp::DiscreteDistribution& rate_dist,
                    bpp::TreeTemplate<bpp::Node>& tree);

    /// Calculate the sum of log-likelihoods. Alias for #log_likelihood
    double operator()();

    /// Calculate the sum of log likelihoods
    double log_likelihood();

    /// Gets the BEAGLE likelihood calculator
    inline std::shared_ptr<Beagle_tree_likelihood> calculator() { return calculator_; }

private:
    std::shared_ptr<Beagle_tree_likelihood> calculator_;
    std::vector<Tree_log_likelihood> additional_log_likes;

    bpp::TreeTemplate<bpp::Node>* tree;
};

}} // Namespace

#endif

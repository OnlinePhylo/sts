#ifndef STS_ONLINE_EXPONENTIAL_BRANCH_LENGTH_PRIOR_H
#define STS_ONLINE_EXPONENTIAL_BRANCH_LENGTH_PRIOR_H

#include <Bpp/Phyl/TreeTemplate.h>
#include <functional>

namespace sts { namespace online {

/// \brief Generic prior density on branch lengths
class Branch_length_prior
{
public:
    /// \brief Construct with a prior density function
    ///
    /// \param log_prior_density a function of one parameter, which given a branch length, returns the prior density.
    Branch_length_prior(std::function<double(double)> log_prior_density);

    double operator()(const bpp::TreeTemplate<bpp::Node>& tree) const;
private:
    std::function<double(double)> log_prior_density;
};

}}  // Namespaces

#endif
